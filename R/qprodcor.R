#' @rdname prodcor
#' @order 3
#' @export
qprodcor <- function( p, rho ) {
    # check inputs
    if ( missing( p ) )
        stop( '`p` is mandatory!' )
    if ( missing( rho ) )
        stop( '`rho` is mandatory!' )
    if ( !is.numeric( p ) )
        stop( '`p` must be numeric!' )
    if ( !is.numeric( rho ) )
        stop( '`rho` must be numeric!' )
    # p can be a vector, but rho cannot
    if ( length( rho ) != 1 )
        stop( '`rho` must be scalar!' )
    # in this function, rho has to be in [-1,1]
    if ( rho > 1 )
        stop( '`rho` must be equal to or less than 1!' )
    if ( rho < -1 )
        stop( '`rho` must be equal to or greater than -1!' )
    # and p has to be a probability
    if ( min( p ) < 0 )
        stop( '`p` must be non-negative!' )
    if ( max( p ) > 1 )
        stop( '`p` must be smaller or equal to 1!' )

    # use chi squared if rho=1
    if ( rho == 1 ) {
        return( stats::qchisq( p, df = 1 ) )
    } else if ( rho == -1 )
        # this is the adequate transformation for negative correlations
        return( - stats::qchisq( 1 - p, df = 1 ) )
    # else continue...
    
    # start calculating individual values
    n <- length( p )
    y <- rep.int( NA, n )
    # this is a coarse approximation useful as a starting point for the root finder, which is particularly tuned for median finding
    # guess quantile under chi-squared (which corresponds to rho=1)
    # reverse quantiles for negative correlations
    # dampen by given rho, also works well when negative
    y2 <- rho * stats::qchisq( if ( rho > 0 ) p else 1 - p, df = 1 )
    
    for ( i in 1 : n ) {
        pi <- p[ i ]

        # handle edge cases
        if ( pi == 1 ) {
            y[ i ] <- Inf
        } else if ( pi == 0 ) {
            y[ i ] <- -Inf
        } else {
            # these bounds work well for median, but in other cases may be off, so allow for interval to be extended (unfortunately you can't just set infinities)
            # add tolerance here too to ensure "opposite signs" at lower and upper
            if ( rho > 0 ) {
                lower <- y2[i] - 0.069
                upper <- y2[i] + 1e-3
            } else {
                lower <- y2[i] - 1e-3
                upper <- y2[i] + 0.069
            }
            y[ i ] <- stats::uniroot( function( z ) pprodcor( z, rho ) - pi, c( lower, upper ), extendInt = 'yes' )$root
        }
    }
    y
}
