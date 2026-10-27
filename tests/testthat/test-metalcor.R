test_that( "This bug is fixed: 'object length is not a multiple of subscript length'", {
    # this very specific, very simple toy example caused a bug in the latest R-devel, before I fixed the bug (now it works fine)
    expect_silent(
        dprodcor( c( -Inf, -1000, -1 ), -0.6 )
    )
})

validate_prodcor <- function( n, x, p, rho ) {
    ### dprodcor

    # create failure on purpose
    expect_error( dprodcor() )
    expect_error( dprodcor( x ) )
    expect_error( dprodcor( rho = rho ) )
    if ( n > 1 )
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
    if ( n > 1 )
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
    if ( n > 1 )
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
    if ( n > 1 )
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

    # a specific relevant edge case, for median calculations
    # backtracking from quantiles to cumulatives
    # these should decrease as rho increases
    y99 <- pprodcor( em, 0.99 )
    y999 <- pprodcor( em, 0.999 )
    y1 <- pprodcor( em, 1 )
    expect_true( y99 > y999 )
    expect_true( y999 > y1 )

    # downstream effect 
    y99 <- qprodcor( 0.5, 0.99 )
    y999 <- qprodcor( 0.5, 0.999 )
    y1 <- qprodcor( 0.5, 1 )
    expect_true( y99 < y999 )
    expect_true( y999 < y1 )
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

    # run a new example from vignette that failed (until debugged)
    x <- ( (-1500) : 1500 ) / 1000
    n <- length( x )
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

validate_metalcor <- function( out, S ) {
    # basic validation of this output
    expect_true( is.list( out ) )
    expect_equal( names( out ), c('assoc', 'R') )
    # validate assoc
    assoc <- out$assoc
    expect_true( is.data.frame( assoc ) )
    expect_equal( colnames( assoc ), c('id', 'chr', 'pos', 'ref', 'alt', 'n', 'n_studies', 'beta', 'se', 'z', 'p') )
    expect_true( all( assoc$n > 0 ) )
    expect_true( all( assoc$n_studies > 0 ) )
    expect_true( all( assoc$n_studies <= S ) )
    expect_true( all( assoc$se > 0 ) )
    expect_true( all( assoc$p >= 0 ) )
    expect_true( all( assoc$p <= 1 ) )
    
    # confirm that data is sorted as desired!
    # this assumes numeric chromosomes!
    expect_true( all( diff( assoc$chr ) >= 0 ) )
    # positions are more complicated, since each chr as-is resets pos to 0. this adds an offset that makes it work
    # first get max position per chromosome
    chr2posmax <- dplyr::summarize( dplyr::group_by( assoc, chr ), posmax = max( pos ) )
    nchr <- nrow( chr2posmax )
    if ( nchr > 1 ) {
        # then get cumulative sum
        chr2posmax$cumposmax <- cumsum( chr2posmax$posmax )
        # want first chromosome to be zero, and the rest to be shifted
        chr2posmax$cumposmax[ 2 : nchr ] <- chr2posmax$cumposmax[ 1 : ( nchr - 1 ) ]
        chr2posmax$cumposmax[ 1 ] <- 0
        assoc$pos2 <- assoc$pos + chr2posmax$cumposmax[ match( assoc$chr, chr2posmax$chr ) ]
    } else
        # single chromosome has trivial shift
        assoc$pos2 <- assoc$pos
    expect_true( all( diff( assoc$pos2 ) >= 0 ) )
    
    # validate R
    R <- out$R
    expect_true( is.matrix( R ) )
    expect_equal( dim( R ), c(S, S) )
    # as a covariance matrix, diagonal is always non-negative
    expect_true( all( diag( R ) > 0 ) )
}

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
        # have more than one chr, but to keep high overlap let's have it be deterministic with pos
        # ref and alt just have to be present, and different...
        ids <- sample.int( n, n_sub )
        study <- data.frame(
            id = ids,
            chr = 1 + round( ids / 250 ),
            pos = ids,
            ref = 'A',
            alt = 'G',
            n = i * 1000,
            beta = rnorm( n_sub ),
            se = abs( rnorm( n_sub ) ) + 1e-8 # should be positive, avoid perfect zeroes
        )
        studies <- c( studies, list( study ) )
    }

    # alternate version with random loci flipped in the 2nd and 3rd studies
    # this will safely edit a copy
    studies2 <- studies
    # don't flip IDs missing in first study, that will change output and complicate comparison
    study1_ids <- studies2[[ 1 ]]$id
    for ( i in 2 : S ) {
        # working copy
        study <- studies2[[ i ]]
        # flip half of each dataset
        indexes <- sample.int( n_sub, n_sub / 2 )
        # keep only loci present in the first study
        indexes <- indexes[ study$id[ indexes ] %in% study1_ids ]
        # change sign of coefficient
        study$beta[ indexes ] <- - study$beta[ indexes ]
        # and switch ref and alt
        tmp <- study$ref[ indexes ]
        study$ref[ indexes ] <- study$alt[ indexes ]
        study$alt[ indexes ] <- tmp
        # put update back in list
        studies2[[ i ]] <- study
    }
    # versions of both inputs without IDs
    studies_noid <- studies
    studies2_noid <- studies2
    for ( i in 1 : S ) {
        studies_noid[[ i ]]$id <- NULL
        studies2_noid[[ i ]]$id <- NULL
    }
    
    for ( median in c(TRUE, FALSE) ) {
        expect_silent(
            out <- metalcor( studies, median = median, useid = TRUE )
        )
        validate_metalcor( out, S )
        
        # repeat for flipped data, should get the same answer!
        expect_silent(
            out2 <- metalcor( studies2, median = median, useid = TRUE )
        )
        expect_equal( out2, out )

        # should also get the same results in these cases when ignoring IDs
        # except output also doesn't have IDs
        out_noid <- out
        out_noid$assoc$id <- NULL
        expect_silent(
            out2 <- metalcor( studies, median = median, useid = FALSE )
        )
        expect_equal( out2, out_noid )
        expect_silent(
            out2 <- metalcor( studies_noid, median = median, useid = FALSE )
        )
        expect_equal( out2, out_noid )

        expect_silent(
            out2 <- metalcor( studies2, median = median, useid = FALSE )
        )
        expect_equal( out2, out_noid )
        expect_silent(
            out2 <- metalcor( studies2_noid, median = median, useid = FALSE )
        )
        expect_equal( out2, out_noid )

        # noswitch also gives same result on aligned `studies`, but they will differ for `studies2`
        expect_silent(
            out2 <- metalcor( studies, median = median, noswitch = TRUE, useid = TRUE )
        )
        expect_equal( out2, out )
        expect_silent(
            out2 <- metalcor( studies, median = median, noswitch = TRUE, useid = FALSE )
        )
        expect_equal( out2, out_noid )
        expect_silent(
            out2 <- metalcor( studies_noid, median = median, noswitch = TRUE, useid = FALSE )
        )
        expect_equal( out2, out_noid )
        
        expect_silent(
            out3 <- metalcor( studies2, median = median, noswitch = TRUE, useid = TRUE )
        )
        validate_metalcor( out3, S )
        out3_noid <- out3
        out3_noid$assoc$id <- NULL
        expect_silent(
            out2 <- metalcor( studies2, median = median, noswitch = TRUE, useid = FALSE )
        )
        expect_equal( out2, out3_noid )
        expect_silent(
            out2 <- metalcor( studies2_noid, median = median, noswitch = TRUE, useid = FALSE )
        )
        expect_equal( out2, out3_noid )
    }
})

