#' @rdname prodcor
#' @order 2
#' @export
pprodcor <- function( x, rho, eps = 1e-8 ) {
    # check inputs
    if ( missing( x ) )
        stop( '`x` is mandatory!' )
    if ( missing( rho ) )
        stop( '`rho` is mandatory!' )
    if ( !is.numeric( x ) )
        stop( '`x` must be numeric!' )
    if ( !is.numeric( rho ) )
        stop( '`rho` must be numeric!' )
    # x can be a vector, but rho cannot
    if ( length( rho ) != 1 )
        stop( '`rho` must be scalar!' )
    # in this function, rho has to be in [-1,1]
    if ( rho > 1 )
        stop( '`rho` must be equal to or less than 1!' )
    if ( rho < -1 )
        stop( '`rho` must be equal to or greater than -1!' )
    # x can be any real number
    if ( eps <= 0 )
        stop( '`eps` must be strictly positive!' )
    
    # use chi squared if rho=1
    if ( rho == 1 ) {
        return( stats::pchisq( x, df = 1 ) )
    } else if ( rho == -1 )
        # this is the adequate transformation for negative correlations
        return( 1 - stats::pchisq( -x, df = 1 ) )
    # else continue...
    
    # start calculating individual values
    n <- length( x )
    y <- rep.int( NA, n )
    for ( i in 1 : n ) {
        xi <- x[ i ]
        if ( xi < -eps ) {
            try( y[ i ] <- stats::integrate( dprodcor, lower = -Inf, upper = xi, rho = rho )$value, silent = TRUE )
            while ( is.na( y[ i ] ) ) {
                # failures are expected near zero, also make bigger until it succeeds
                xi <- xi * 10
                try( y[ i ] <- stats::integrate( dprodcor, lower = -Inf, upper = xi, rho = rho )$value, silent = TRUE )
            }
        } else if ( xi < eps ) {
            # this case ( xi ~ 0 ) fails often when rho is large :(
            eps2 <- eps
            try( y[ i ] <- stats::integrate( dprodcor, lower = -Inf, upper = -eps2, rho = rho )$value, silent = TRUE )
            while ( is.na( y[ i ] ) ) {
                # make eps bigger until we get a good number
                eps2 <- eps2 * 10
                try( y[ i ] <- stats::integrate( dprodcor, lower = -Inf, upper = -eps2, rho = rho )$value, silent = TRUE )
            }
        } else {
            try( y[ i ] <- 1 - stats::integrate( dprodcor, lower = xi, upper = Inf, rho = rho )$value, silent = TRUE )
            while ( is.na( y[ i ] ) ) {
                # failures are expected near zero, also make bigger until it succeeds
                xi <- xi * 10
                try( y[ i ] <- 1 - stats::integrate( dprodcor, lower = xi, upper = Inf, rho = rho )$value, silent = TRUE )
            }
        }
    }
    y
}
