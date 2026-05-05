# constant, marks transition from covariance to variance
em <- qchisq( 0.5, 1 ) # [1] 0.4549364

#' Estimated Correlations from Sample Medians of the Correlated Standard Normal Product Distribution
#'
#' This code implements estimation from the sample median of the correlation of a product of two correlated standard normal variables.
#' While this estimation is straightforward from the mean, using the median requires the transformation implemented here because the distribution has heavy tails, so for positive correlations the median is much smaller than the mean.
#' The median version is desired due to its robustness to outliers.
#'
#' The highest correlation of `rho = 1` results in a median of only `q~0.455`.  To permit application to arbitrary data, including data that does not satisfy assumptions and is therefore inflated, the function returns `x/q` for elements whose absolute values of `x/q` exceed 1, and can therefore return values that exceed \[-1, 1\].
#'
#' @param x vector of sample medians
#'
#' @return Vector of estimated correlations, extended for inflated cases.
#'
#' @examples
#'
#' # a large number of samples is required for the median to yield an accurate estimate
#' n <- 10000
#' rho <- 0.3
#' x <- median( rprodcor( n, rho ) )
#' rho_from_median( x )
#'
#' @seealso [qprodcor()]
#' 
#' @export
rho_from_median <- function( x ) {
    # check inputs
    if ( missing( x ) )
        stop( '`x` is mandatory!' )
    if ( !is.numeric( x ) )
        stop( '`x` must be numeric!' )
    # x can otherwise be a vector with any real values...
    
    n <- length( x )
    y <- rep.int( NA, n )
    for ( i in 1 : n ) {
        # solutions are close to this in a sense, and plot shows it as a useful bound (sign dependent)
        rho0 <- x[ i ] / em
        if ( abs( rho0 ) >= 1 - 1e-3 ) {
            # desired value for variance case, which is only case expected to give huge values
            y[ i ] <- rho0
        } else {
            # add tolerance here too to ensure "opposite signs" at lower and upper
            if ( rho0 > 0 ) {
                lower <- rho0 - 1e-3
                upper <- min( 1, rho0 + 0.15 )
            } else {
                lower <- max( -1, rho0 - 0.15 )
                upper <- rho0 + 1e-3
            }
            y[ i ] <- stats::uniroot( function( rho ) qprodcor( 0.5, rho ) - x[ i ], c(lower, upper) )$root
        }
    }
    y
}
