#' Estimate Covariance Matrix of Study Z-scores
#'
#' The default "median" estimator is expected to be more robust to outliers (present due to causal variants and linkage disequilibrium).
#' For completeness, the more standard "mean" estimator is also implemented, which is better when there are no outliers.
#' The determinant is also tested, and if its absolute value is less than `tol` an error is triggered, warning the user that perhaps the same study was included multiple times (as this is otherwise expected to be extremely rare).
#'
#' @param Z Matrix of z-scores, loci along rows, `S` studies along columns.  Missing values can be present and are appropriately ignored.
#' @param median If `TRUE`, uses [rho_from_median()] to calculate covariance values from sample medians; otherwise, sample means are used.
#' @param tol Determinant tolerance for singularity test
#'
#' @return The `S`-by-`S` covariance matrix.
#'
#' @examples
#'
#' # simulate some correlated Z-scores for 3 studies and `m` loci
#' m <- 10000
#' # true correlation matrix
#' R <- matrix(
#'   c( 1.0, 0.2, 0.0,
#'      0.2, 1.0, 0.1,
#'      0.0, 0.1, 1.0 ),
#'   nrow = 3, ncol = 3
#' )
#' # use Cholesky decomposition to simulate via affine transform
#' L <- chol( R )
#' Z <- matrix( rnorm( m * 3 ), nrow = m, ncol = 3 ) %*% L
#'
#' # calulate estimates
#' estimate_R( Z )
#' estimate_R( Z, median = FALSE )
#' # compare to naive covariance estimate, which does not assume zero mean for z-scores
#' cov( Z, use = 'pairwise.complete.obs' )
#'
#' @seealso [rho_from_median()]
#' 
#' @export
estimate_R <- function ( Z, median = TRUE, tol = 1e-8 ) {
    # check inputs
    if ( missing( Z ) )
        stop( '`Z` is mandatory!' )
    if ( !is.matrix( Z ) )
        stop( '`Z` must be a matrix!' )
    if ( !is.numeric( Z ) )
        stop( '`Z` must be numeric!' )
    
    # estimate z-score correlation matrix
    if ( median ) {
        # unfortunately we have to build it more carefully...
        S <- ncol( Z ) # number of studies
        R <- matrix( NA, S, S )
        for ( j in 1 : S ) {
            zj <- Z[ , j ]
            for ( k in 1 : j ) {
                mjk <- median( zj * Z[ , k ], na.rm = TRUE )
                rjk <- rho_from_median( mjk )
                R[ j, k ] <- rjk
                R[ k, j ] <- rjk
            }
        }
    } else {
        # this forms the correlation estimate, handles missingness but overfits means
        #R <- cov( Z, use = 'pairwise.complete.obs' )
        
        # this version assumes null more properly, but does not handle missingness
        #R <- crossprod( Z ) / nrow( Z )
        
        # when there's missingness, this is more painful
        # first, number of complete observations
        M <- crossprod( !is.na( Z ) )
        # now missing cases can be zeroed
        Z[ is.na( Z ) ] <- 0
        # we can now compute the final estimate
        R <- crossprod( Z ) / M
    }

    # R is going to have to be invertible... check determinant now
    # (sample estimates on reasonable data are extremely unlikely to cause problems, data essentially has to be duplicated for problems to occur)
    if ( abs( det( R ) ) < tol )
        stop( 'Z-score correlation matrix is numerically singular!  Did you input the same study more than one time?  Det = ', det( R ) )
    
    R
}
