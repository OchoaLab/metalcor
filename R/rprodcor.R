#' @rdname prodcor
#' @order 4
#' @export
rprodcor <- function( n, rho ) {
    # check inputs
    if ( missing( n ) )
        stop( '`n` is mandatory!' )
    if ( missing( rho ) )
        stop( '`rho` is mandatory!' )
    if ( !is.numeric( n ) )
        stop( '`n` must be numeric!' )
    if ( !is.numeric( rho ) )
        stop( '`rho` must be numeric!' )
    # if n is a vector, replace it with its length, as rnorm does
    if ( length( n ) > 1 ) {
        n <- length( n )
    } else {
        # otherwise make sure it's a positive integer
        n <- as.integer( n )
        if ( n < 1 )
            stop( '`n` must round to a positive value!' )
    }
    # x can be a vector, but rho cannot
    if ( length( rho ) != 1 )
        stop( '`rho` must be scalar!' )
    # in this function, rho has to be in [-1,1]
    if ( rho > 1 )
        stop( '`rho` must be equal to or less than 1!' )
    if ( rho < -1 )
        stop( '`rho` must be equal to or greater than -1!' )

    # use chi squared if rho=1
    if ( rho == 1 ) {
        return( stats::rchisq( n, df = 1 ) )
    } else if ( rho == -1 )
        # this is the adequate transformation for negative correlations
        return( - stats::rchisq( n, df = 1 ) )
    # else continue...
    
    # the goal is to generate two random standard normal variables that are correlated, then multiply them
    # use affine transformation of the covariance matrix R
    # R <- matrix( c( 1, rho, rho, 1 ), 2, 2)
    # cholesky has closed form! (this is transposed from what `chol` returns, but it is what we want for this transform)
    U <- matrix( c( 1, rho, 0, sqrt(1 - rho^2) ), 2, 2 )
    
    # draw the pair of independent standard normals, arranged as a wide matrix with two rows
    # multiplying by U matrix gives us the desired correlated variables
    X <- U %*% matrix( stats::rnorm( 2 * n ), nrow = 2, ncol = n )
    # lastly, multiply the two values of each column and return, as a regular vector
    drop( X[1,] * X[2,] )
}
