validate_prodcor <- function( n, x, p, rho ) {
    ### dprodcor

    # create failure on purpose
    expect_error( dprodcor() )
    expect_error( dprodcor( x ) )
    expect_error( dprodcor( rho = rho ) )
    expect_error( dprodcor( x, x ) )

    # a successful case
    expect_silent( 
        y <- dprodcor( x, rho )
    )
    expect_true( is.numeric( y ) )
    expect_true( !anyNA( y ) )
    expect_equal( length( y ), n )
    expect_true( min( y ) >= 0 )

    ### pprodcor

    # create failure on purpose
    expect_error( pprodcor() )
    expect_error( pprodcor( x ) )
    expect_error( pprodcor( rho = rho ) )
    expect_error( pprodcor( x, x ) )
    
    # a successful case
    expect_silent( 
        y <- pprodcor( x, rho )
    )
    expect_true( is.numeric( y ) )
    expect_true( !anyNA( y ) )
    expect_equal( length( y ), n )
    expect_true( min( y ) >= 0 )
    # here additionally we expect values to increase monotonically, but we need to order by x
    expect_true( all( diff( y[ order( x ) ] ) >= 0 ) )

    ### qprodcor

    # create failure on purpose
    expect_error( qprodcor() )
    expect_error( qprodcor( p ) )
    expect_error( qprodcor( rho = rho ) )
    expect_error( qprodcor( p, p ) )
    
    # a successful case
    expect_silent( 
        y <- qprodcor( p, rho )
    )
    expect_true( is.numeric( y ) )
    expect_true( !anyNA( y ) )
    expect_equal( length( y ), n )
    # quantiles have unlimited range

    # also test p of 0 or 1 specifically
    expect_silent( 
        y <- qprodcor( p = c(0,1), rho )
    )
    # limits are all reals, except in two edge cases
    y_exp <- c( if ( rho == 1 ) 0 else -Inf, if ( rho == -1 ) 0 else Inf)
    expect_equal( y, y_exp )

    ### rprodcor

    # create failure on purpose
    expect_error( rprodcor() )
    expect_error( rprodcor( n ) )
    expect_error( rprodcor( rho = rho ) )
    expect_error( rprodcor( n, n ) )
    
    # a successful case
    expect_silent( 
        y <- rprodcor( n, rho )
    )
    expect_true( is.numeric( y ) )
    expect_true( !anyNA( y ) )
    expect_equal( length( y ), n )
    # quantiles have unlimited range
}

test_that("*prodcor family works", {
    # simulate some data
    n <- 100
    x <- rnorm( n )
    p <- runif( n )

    # validate all cases for a random rho
    validate_prodcor( n, x, p, rho = runif( 1, -1, 1 ) )
    # some special/edge cases
    validate_prodcor( n, x, p, rho = 1 )
    validate_prodcor( n, x, p, rho = -1 )
    # not really a "edge" case, but still expected to be a very common request
    validate_prodcor( n, x, p, rho = 0 )

})

test_that( "rho_from_median works", {
    # size of window around troublesome edge cases
    eps <- 1e-4
    # create these many samples from each case (4 cases), plus 3 fixed cases
    nc <- 10

    # these are all the cases to test!
    x <- c(
        # fixed cases
        -1, 0, 1,
        # test random values near 1
        em * runif( nc, 1 - eps, 1 + eps ),
        # test random values near -1
        - em * runif( nc, 1 - eps, 1 + eps ),
        # test random values near 0
        runif( nc, -eps, eps ),
        # test totally random numbers
        runif( nc, -1, 1 )
    )
    # actual length
    n <- length( x )

    # cause errors on purpose
    expect_error( rho_from_median( ) )
    expect_error( rho_from_median( letters ) )

    # a succesful case
    expect_silent( 
        y <- rho_from_median( x )
    )
    expect_true( is.numeric( y ) )
    expect_true( !anyNA( y ) )
    expect_equal( length( y ), n )
    # these should have the same sign
    expect_true( min( y * x ) >= 0 )
})

test_that( "estimate_R works", {
    # dimensions of simulation
    m <- 10000
    # true correlation matrix
    R <- matrix(
        c(
            1.0, 0.2, 0.0,
            0.2, 1.0, 0.1,
            0.0, 0.1, 1.0
        ),
        nrow = 3, ncol = 3
    )
    # use Cholesky decomposition to simulate via affine transform
    L <- chol( R )
    Z <- matrix( rnorm( m * 3 ), nrow = m, ncol = 3 ) %*% L

    # calulate estimates
    for ( median in c(TRUE, FALSE) ) {
        expect_silent(
            R_est <- estimate_R( Z, median = median )
        )
        expect_true( is.matrix( R_est ) )
        expect_true( is.numeric( R_est ) )
        expect_equal( dim( R_est ), c(3, 3) )
        expect_true( min( diag( R_est ) ) >= 0 )
    }
})



test_that( "metalcor works", {
    # pick a dense subset of 1000 loci, with some non-overlap to test missingness properly
    n <- 1000
    n_sub <- 950
    # number of studies
    S <- 3

    # construct studies
    studies <- NULL
    for ( i in 1 : S ) {
        # ids have to be unique
        # have more than one chr, but to keep high overlap let's hav eit be deterministic with pos
        ids <- sample.int( n, n_sub )
        study <- data.frame(
            id = ids,
            chr = 1 + round( ids / 250 ),
            pos = ids,
            n = i * 1000,
            beta = rnorm( n_sub ),
            se = rnorm( n_sub )
        )
        studies <- c( studies, list( study ) )
    }

    for ( median in c(TRUE, FALSE) ) {
        expect_silent(
            out <- metalcor( studies, median = median )
        )
    }
})

