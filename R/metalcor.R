#' Meta-Analysis of Correlated Genetic Association Studies
#'
#' This function accepts the summary statistics from several studies, estimates their z-score correlation matrix, and performs the meta-analysis taking those correlations into account.
#'
#' @param studies A list of `S` studies, each of which is a data frame with the following required columns: `id`, `chr`, `pos`, `beta`, `se`, `n`.
#' Column names have to match exactly.
#' `id` is used to match variants between studies (they must be unique for different variants, and be the same for the same variant).
#' `chr` and `pos` are used solely to sort variants after the studies have been merged.
#' `beta` and `se` are required to perform meta-analysis.
#' `n` is used to report locus-specific total sample sizes under missingness.
#' @param median If `TRUE` (default), estimates correlations from sample medians instead of means.
#' Passed to [estimate_R()].
#' @param df Degree of freedom for calculating the p-values from the z-scores (default of 1 is best practically always).
#' @param tol Determinant tolerance for singularity test of R estimate
#'
#' @return A list with two named elements:
#' - `assoc`: the meta-analyzed study, a tibble with columns `id`, `chr`, `pos`, `n`, `n_studies`, `beta`, `se`, `z`, and `p`; sorted by `chr` and `pos`.
#' - `R`: the `S`-by-`S` estimated z-score covariance matrix
#'
#' @examples
#' # construct two toy studies just to run example, with minimal columns
#' study1 <- data.frame(
#'   id = paste0( 'rs', 1:5 ),
#'   chr = 1,
#'   pos = 1:5,
#'   n = 2000,
#'   beta = rnorm( 5 ),
#'   se = rnorm( 5 )
#' )
#' # note the second study is missing the 5th SNP, this is fine
#' study2 <- data.frame(
#'   id = paste0( 'rs', 1:4 ),
#'   chr = 1,
#'   pos = 1:4,
#'   n = 5000,
#'   beta = rnorm( 4 ),
#'   se = rnorm( 4 )
#' )
#'
#' # gather the studies in a list
#' studies <- list( study1, study2 )
#' # this performs the meta-analysis modeling covariance!
#' out <- metalcor( studies )
#' # this is the meta-analyzed association table
#' out$assoc
#' # and this is the estimated study covariance matrix
#' out$R
#'
#' @seealso [estimate_R()]
#' 
#' @export
metalcor <- function( studies, median = TRUE, df = 1, tol = 1e-8 ) {
    # check inputs
    if ( missing( studies ) )
        stop( '`studies` is mandatory!' )
    if ( !is.list( studies ) )
        stop( '`studies` must be a list!' )
    S <- length( studies )
    if ( S < 2 )
        stop( 'Provide at least two studies!' )
    
    # this merges the studies!
    study_full <- merge_studies( studies )

    # pull out the desired matrices to pass to the math code
    betas <- as.matrix( dplyr::select( study_full, tidyselect::num_range( "beta", 1:S ) ) )
    ses <- as.matrix( dplyr::select( study_full, tidyselect::num_range( "se", 1:S ) ) )
    # colnames don't add anything here, better to remove them for clarity
    colnames( betas ) <- NULL
    colnames( ses ) <- NULL
    # exclude the separate betas and ses from the data frames now
    study_full <- dplyr::select( study_full, 'id', 'chr', 'pos', 'n' )

    # estimate the correlation matrix
    # compute z scores matrix on the fly, by elementwise normalization
    R <- estimate_R( betas / ses, median = median, tol = tol )
    
    # this performs the actual calculations
    gwas_meta <- meta_corr_math( betas, ses, R, df = df )

    # add back the id/chr/pos and we're done
    gwas_meta <- dplyr::bind_cols( study_full, gwas_meta )
    # return R matrix too
    list(
        assoc = gwas_meta,
        R = R
    )
}

merge_studies <- function( studies ) {
    S <- length( studies )

    # check that each element is a data frame, errors get confusing otherwisee
    for ( j in 1 : S )
        if ( !is.data.frame( studies[[ j ]] ) )
            stop( 'Study ', j, ' was not a data frame!' )
        
    # make a merged data frame, doing full join of two studies at the time
    # 1) each study we toss the non-essential columns,
    # 2) columns to merge on (id, chr, pos) are kept same name in each study pre-merge, but columns to keep different per study are explicitly numbered
    study_full <- dplyr::select( studies[[1]], 'id', 'chr', 'pos', beta1 = 'beta', se1 = 'se', n1 = 'n' )
    for ( j in 2 : S ) {
        # use only desired columns
        study_j <- dplyr::select( studies[[j]], 'id', 'chr', 'pos', 'beta', 'se', 'n' )
        # append number to columns we want differentiated
        study_j <- dplyr::rename_with( study_j, ~ paste0(.x, j), .cols = c( 'beta', 'se', 'n' ) )
        # now we can apply the join
        study_full <- dplyr::full_join( study_full, study_j, by = dplyr::join_by( 'id', 'chr', 'pos' ) )
    }

    # trivial sum of sample sizes across non-missing studies, no other use for the n's
    study_full$n <- rowSums( dplyr::select( study_full, tidyselect::num_range( "n", 1:S ) ), na.rm = TRUE )
    # immediately remove the separate n values
    study_full <- dplyr::select( study_full, !tidyselect::num_range( "n", 1:S ) )

    # data can be kinda sorted but not quite, lets sort now (only use for chr/pos) 
    study_full <- dplyr::arrange( study_full, 'chr', 'pos' )

    # done merging!
    study_full
}

# this does the actual math, separated from the more tedious alignment code
meta_corr_math <- function( betas, ses, R, df = 1 ) {
    # useful to know the total number of studies
    S <- ncol( betas )
    # count number of studies with observations per locus, nice stat and helps a bit with code
    Ss <- rowSums( !is.na( betas ) )

    # construct the sigma of every locus explicitly, appears to behave way better numerically, but should be slower
    m <- nrow( betas )
    ses_meta <- rep.int( NA, m )
    betas_meta <- rep.int( NA, m )
    for ( i in 1 : m ) {
        # deal with missingness
        # we pre-calculated numbers of observations per locus, pull that now here
        S_obs <- Ss[ i ]
        # can't think of a reasonable scenario where all data is missing, but just handle it anyway
        # the output vectors are filled with NAs, if we don't touch them we're good in this case
        if ( S_obs == 0 )
            next
        # pull vectors now, we'll need them
        ses_i <- ses[ i, ]
        betas_i <- betas[ i, ]
        R_i <- R # a copy, in case it needs subsetting
        if ( S_obs < S ) {
            # subset everything now, since it's needed
            # se, beta will be missing together under any reasonable scenario
            indexes <- !is.na( ses_i )
            ses_i <- ses_i[ indexes ]
            betas_i <- betas_i[ indexes ]
            # shortcut analysis if there's just one non-missing study
            if ( S_obs == 1 ) {
                # just copy the one non-missing value
                ses_meta[ i ] <- ses_i
                betas_meta[ i ] <- betas_i
                # move on to the next SNP! 
                next
            }
            # only in this case it's worth subsetting
            R_i <- R_i[ indexes, indexes ]
        }
        # form explicit beta covariance matrix, using correlation calculated earlier
        D <- diag( ses_i )
        # try to invert, though in some cases this fails.  This is always because some ses_i values are very close to zero (R_i is otherwise well behaved in practice, and was checked earlier anyway)
        Sigma_inv <- try( solve( D %*% R_i %*% D ), silent = TRUE )
        if ( "try-error" %in% class( Sigma_inv ) ) {
            # again, one sigma is extremely close to zero...
            # first let's identify it
            j <- which.min( ses_i )
            # in these cases precision is kinda pointless, asymptotics show the beta and se of that worse case is roughly the value approached (with some annoying constants we just ignore here)
            # if beta=0 but se>0, below we'll calculate a p-value of 1 (this is perfect)
            # if se=0 and beta>0 (I think this unlikely), we'd get z=Inf and p=0, let's avoid that!!!  Also cap extreme z-scores (bigger than 4 in this case), force beta to zero in those cases so z score is also zero.
            # lastly, se=0 and beta=0 gives z=NaN and p=NaN too, let's keep that too and hope it's not too unpleasant downstream (at least for my evals (aucpr, srmsdp, infl) this is fine)
            ses_meta[ i ] <- ses_i[ j ]
            betas_meta[ i ] <- if ( ses_i[ j ] == 0 || betas_i[ j ] / ses_i[ j ] > 4 ) 0 else betas_i[ j ]
        } else {
            # regular handling, most common
            denom_i <- sum( Sigma_inv )
            ses_meta[ i ] <- 1 / sqrt( denom_i )
            betas_meta[ i ] <- sum( Sigma_inv %*% betas_i ) / denom_i
        }
    }
    
    # meta z-scores!
    zs_meta <- betas_meta / ses_meta
    # p-values
    ps_meta <- stats::pchisq( zs_meta^2, df, lower.tail = FALSE )
    # organize into a tibble and return
    gwas_meta <- tibble::tibble(
                             n_studies = Ss,
                             beta = betas_meta,
                             se = ses_meta,
                             z = zs_meta,
                             p = ps_meta
                         )
}
