#' Meta-Analysis of Correlated Genetic Association Studies
#'
#' This function accepts the summary statistics from several studies, estimates their z-score correlation matrix, and performs the meta-analysis taking those correlations into account.
#'
#' @param studies A list of `S` studies, each of which is a data frame with the following required columns: `chr`, `pos`, `ref`, `alt`, `beta`, `se`, `n`.
#' Column names have to match exactly.
#' `chr`, `pos`, `ref` and `alt` are jointly used to match variants between studies (they must be unique for different variants, and be the same for the same variant), except if `noswitch = FALSE` (default) then `ref` and `alt` values may be switched in different studies.
#' `chr` and `pos` are also used to sort variants after the studies have been merged.
#' `beta` and `se` are required to perform meta-analysis.
#' `n` is used to report locus-specific total sample sizes under missingness.
#' @param median If `TRUE` (default), estimates correlations from sample medians instead of means.
#' Passed to [estimate_R()].
#' @param df Degree of freedom for calculating the p-values from the z-scores (default of 1 is best practically always).
#' @param tol Determinant tolerance for singularity test of R estimate
#' @param useid Additionally use IDs for matching and in outputs.
#' If `FALSE` (default), IDs are totally ignored, even if inputs had them.
#' If `TRUE`, `id` will be used to match variants (together with `chr`, `pos`, `ref`, `alt`), and there will be `id` in the output.
#' IDs can be useful when ref/alt switches are allowed, to better disambiguate matching allele pairs, especially for unstandardized indels.
#' @param noswitch Do not allow `ref` and `alt` switching.  Default is to attempt to switch them for a locus and study if matching failed otherwise.
#' `noswitch = TRUE` can be a good idea if you're sure most or all of your data is aligned to reference, or when insertions and deletions both appear and are symmetric (for example, in one allele pair "A/AB" and another "AB/A"), which can result in confusion.
#'
#' @return A list with two named elements:
#' - `assoc`: the meta-analyzed study, a tibble with columns `id` (if `useid = TRUE`), `chr`, `pos`, `ref`, `alt`, `n`, `n_studies`, `beta`, `se`, `z`, and `p`; sorted by `chr` and `pos`.
#'   `ref` and `alt` values are those of the first study in which each variant appears, and if a subsequent study has then switched, the variant is aligned by changing the sign of its `beta` prior to meta-analysis.
#' - `R`: the `S`-by-`S` estimated z-score covariance matrix
#'
#' @examples
#' # construct two toy studies just to run example, with minimal columns
#' study1 <- data.frame(
#'   chr = 1,
#'   pos = 1:5,
#'   ref = 'A',
#'   alt = 'G',
#'   n = 2000,
#'   beta = rnorm( 5 ),
#'   se = abs( rnorm( 5 ) )
#' )
#' # note the second study is missing the 5th SNP, this is fine
#' # also, the 4th SNP is has switched ref/alt values
#' study2 <- data.frame(
#'   chr = 1,
#'   pos = 1:4,
#'   ref = c( 'A', 'A', 'A', 'G' ),
#'   alt = c( 'G', 'G', 'G', 'A' ),
#'   n = 5000,
#'   beta = rnorm( 4 ),
#'   se = abs( rnorm( 4 ) )
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
metalcor <- function( studies, median = TRUE, df = 1, tol = 1e-8, useid = FALSE, noswitch = FALSE ) {
    # check inputs
    if ( missing( studies ) )
        stop( '`studies` is mandatory!' )
    if ( !is.list( studies ) )
        stop( '`studies` must be a list!' )
    S <- length( studies )
    if ( S < 2 )
        stop( 'Provide at least two studies!' )
    
    # this merges the studies!
    study_full <- merge_studies( studies, useid, noswitch )

    # pull out the desired matrices to pass to the math code
    betas <- as.matrix( dplyr::select( study_full, tidyselect::num_range( "beta", 1:S ) ) )
    ses <- as.matrix( dplyr::select( study_full, tidyselect::num_range( "se", 1:S ) ) )
    # colnames don't add anything here, better to remove them for clarity
    colnames( betas ) <- NULL
    colnames( ses ) <- NULL
    # exclude the separate betas and ses from the data frames now
    study_full <- dplyr::select( study_full, if (useid) 'id', 'chr', 'pos', 'ref', 'alt', 'n' )

    # estimate the correlation matrix
    # compute z scores matrix on the fly, by elementwise normalization
    R <- estimate_R( betas / ses, median = median, tol = tol )
    
    # this performs the actual calculations
    gwas_meta <- meta_corr_math( betas, ses, R, df = df )

    # add back the id/chr/pos/ref/alt and we're done
    gwas_meta <- dplyr::bind_cols( study_full, gwas_meta )
    # return R matrix too
    list(
        assoc = gwas_meta,
        R = R
    )
}

#' @importFrom rlang .data
merge_studies <- function( studies, useid, noswitch ) {
    S <- length( studies )

    # check that each element is a data frame, errors get confusing otherwise
    for ( j in 1 : S )
        if ( !is.data.frame( studies[[ j ]] ) )
            stop( 'Study ', j, ' was not a data frame!' )

    # if we're not using existing IDs, make temporary IDs that capture chr:pos:ref:alt except ref/alt are sorted if switching is allowed.  IDs are required to facilitate ref/alt switching, to match loci
    if ( !useid && !noswitch )
        for ( j in 1 : S )
            # ID is chr:pos:ref:alt except ref and alt are sorted alphabetically
            studies[[ j ]] <- dplyr::mutate( studies[[ j ]], id = paste( .data$chr, .data$pos, pmin( .data$ref, .data$alt ), pmax( .data$ref, .data$alt ), sep = ":" ) )

    # use ID if requested or when switching
    join_cols <- c( 'chr', 'pos', 'ref', 'alt' )
    if ( useid || !noswitch )
        join_cols <- c( 'id', join_cols )
    
    # make a merged data frame, doing full join of two studies at the time
    # 1) each study we toss the non-essential columns,
    # 2) columns to merge on (id, chr, pos, ref, alt) are kept same name in each study pre-merge, but columns to keep different per study are explicitly numbered
    study_full <- dplyr::select( studies[[1]], tidyselect::all_of( c(join_cols, beta1 = 'beta', se1 = 'se', n1 = 'n') ) )
    for ( j in 2 : S ) {
        # use only desired columns
        study_j <- dplyr::select( studies[[j]], tidyselect::all_of( c(join_cols, 'beta', 'se', 'n') ) )

        # skip all of this if we don't allow switches
        if ( !noswitch ) {
            # identify data that would match except ref/alt are switched
            # first exclude matches by id/chr/pos/ref/alt (we won't alter those), leaving those we might alter
            # NOTE: IDs are always present when we allow switching (constructed if needed)
            tmp <- dplyr::anti_join( study_j, study_full, by = join_cols )
            # in this subset, switch ref/alt now
            tmp <- dplyr::rename( tmp, ref = 'alt', alt = 'ref' )
            # we want the subset that matches now
            tmp <- dplyr::semi_join( tmp, study_full, by = join_cols )
            if ( nrow( tmp ) > 0 ) {
                # now we know which rows to switch, get their indexes
                # here use simple matching by ID
                # if IDs are missing/ignored and we allow switching, we have already constructed IDs that will work
                indexes <- study_j$id %in% tmp$id
                # most important transformation is changing the sign of the coefficient!  No other stats need to change
                study_j$beta[ indexes ] <- - study_j$beta[ indexes ]
                # but we also need to switch ref and alt here, only on these rows
                tmp2 <- study_j$ref[ indexes ]
                study_j$ref[ indexes ] <- study_j$alt[ indexes ]
                study_j$alt[ indexes ] <- tmp2
                # done!
            }
        }
        
        # append number to columns we want differentiated
        study_j <- dplyr::rename_with( study_j, ~ paste0(.x, j), .cols = c( 'beta', 'se', 'n' ) )
        # now we can apply the join, keeping unmatched things as separate rows
        study_full <- dplyr::full_join( study_full, study_j, by = join_cols )
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
