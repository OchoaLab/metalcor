#' The Correlated Standard Normal Product Distribution
#'
#' Density, cumulative, quantile, and random generation for the product of two correlated standard normal variables with correlation parameter `rho`.
#'
#' `pprodcor` does not have a closed form, so it is calculated by numerical integration of `dprodcor`.
#' Furthermore, `dprodcor(0)` is infinite, so in order for numerical integration to succeed, we exclude the window `c(-eps, eps)` for arguments near zero.  `eps` is the initial value, but if failure is still encountered, `eps` is incremented by factors of 10 untill successful.  (For positive arguments integration from above is used instead, subtracted from 1.)
#' `qprodcor` also does not have a closed form, so it is calculated using a root finder on `pprodcor`, which makes it very slow.
#'
#' @param x Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.  If `length(n) > 1`, the length is taken to be the number required.
#' @param rho Correlation parameter, must be a scalar in \[-1, 1\].
#' @param eps Near zero value to use for integrating around the divergence at `x=0`.
#'
#' @return `dprodcor` gives the density, `pprodcor` the cumulative, and `qprodcor` the quantile function.  `rprodcor` generates random deviates.  The length of the result is determined by `n` for `rnorm`, and it is the length of `x` or `p` for the other functions.
#'
#' @references Nadarajah, Saralees, and Tibor K. Pogány. “On the Distribution of the Product of Correlated Normal Random Variables.” Comptes Rendus Mathematique 354, no. 2 (2016): 201–4. https://doi.org/10.1016/j.crma.2015.10.019.
#'
#' Gaunt, Robert E. “The Basic Distributional Theory for the Product of Zero Mean Correlated Normal Random Variables.” Statistica Neerlandica 76, no. 4 (2022): 450–70. https://doi.org/10.1111/stan.12267.
#'
#' @examples
#' n <- 10
#' x <- rnorm( n )
#' p <- runif( n )
#' rho <- 0.8
#' dprodcor( x, rho )
#' pprodcor( x, rho )
#' qprodcor( p, rho )
#' rprodcor( n, rho )
#'
#' @name prodcor
#' @order 1
#' @export
dprodcor <- function( x, rho ) {
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

    # to handle negative rho best in final approx, take advantage of symmetry, assume positive rho onward:
    if ( rho < 0 ) {
        rho <- - rho
        x <- -x
    }
    
    # handle some trivial cases (rho == -1 handled by above symmetry transform)
    if ( rho == 1 )
        return( stats::dchisq( x, df = 1 ) )
    # else continue
    
    # this is the determinant of R
    d <- 1 - rho^2
    # function is in terms of this converted variable
    x <- x / d
    # density from Nadarajah 16
    y <- 1 / ( pi * sqrt( d ) ) * exp( rho * x ) * besselK( abs(x), nu = 0 )
    # just recalculate anything that evaluated to NA with the alternate formula (not sure how worth it it is to find thresholds to decide
    indexes <- is.na( y )
    if ( any( indexes ) ) {
        # handle both positive and negative cases...
        indexes_neg <- x[ indexes ] < 0
        y[ indexes[ indexes_neg ] ] <- 0
        # subset positive cases
        indexes <- indexes[ !indexes_neg ]
        x <- x[ indexes ]
        y[ indexes ] <- 1 / sqrt( 2 * pi * x ) * exp( - x / ( 1 + rho ) )
    }
    # done!
    y
}
