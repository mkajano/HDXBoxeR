## internal helpers for compare_states_unmatched() - not exported

.hl_shift <- function(x, y) {
  ## Hodges-Lehmann shift estimator: median of every pairwise difference
  ## y_j - x_i. The estimand the Wilcoxon rank-sum test is actually a test
  ## of, so it stays consistent with the "classic"/"peptide" p-value story
  ## even when reported alongside a resampled CI.
  stats::median(as.numeric(outer(y, x, "-")))
}

.perm_test_hl <- function(x, y, n_perm) {
  obs <- .hl_shift(x, y)
  pooled <- c(x, y); n1 <- length(x); n2 <- length(y)
  perm_stats <- vapply(seq_len(n_perm), function(i) {
    idx <- sample.int(n1 + n2)
    .hl_shift(pooled[idx[seq_len(n1)]], pooled[idx[(n1 + 1):(n1 + n2)]])
  }, numeric(1))
  (1 + sum(abs(perm_stats) >= abs(obs), na.rm = TRUE)) / (1 + n_perm)
}

## per-peptide matched Welch t-test (replicate-level) for every peptide_id
## present in BOTH s1 and s2 within one region's raw (unaggregated) rows -
## the diagnostic that catches a real but spatially localized effect getting
## diluted to non-significance by the pooled region-level test above: see
## compare_states_unmatched()'s "frac_peptides_significant" section.
.peptide_matched_pvals <- function(bd_region, s1, s2) {
  ids1 <- unique(bd_region$peptide_id[bd_region$Protein.State == s1])
  ids2 <- unique(bd_region$peptide_id[bd_region$Protein.State == s2])
  shared_ids <- intersect(ids1, ids2)
  if (length(shared_ids) == 0) {
    return(data.frame(peptide_id = character(0), p_value = numeric(0), stringsAsFactors = FALSE))
  }
  stats_one <- function(pid) {
    x <- bd_region$deut_value[bd_region$Protein.State == s1 & bd_region$peptide_id == pid]
    y <- bd_region$deut_value[bd_region$Protein.State == s2 & bd_region$peptide_id == pid]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    p <- if (length(x) < 2 || length(y) < 2) NA_real_ else
      tryCatch(stats::t.test(x, y)$p.value, error = function(e) NA_real_)
    c(p_value = if (is.finite(p)) p else NA_real_,
      delta = if (length(x) && length(y)) mean(y) - mean(x) else NA_real_,
      sd1 = if (length(x) >= 2) stats::sd(x) else NA_real_, n1 = length(x),
      sd2 = if (length(y) >= 2) stats::sd(y) else NA_real_, n2 = length(y))
  }
  st <- do.call(rbind, lapply(shared_ids, stats_one))
  data.frame(peptide_id = shared_ids, st, stringsAsFactors = FALSE, row.names = NULL)
}

#' Compare pooled peptide-level HDX uptake across states without matching peptides
#'
#' A deliberately different tool from [compare_states_tp()]: instead of
#' requiring the same peptide to be identified in every state being
#' compared, this pools *every* peptide in a chosen residue range within
#' each state into one distribution, and asks whether that region's typical
#' exchange level differs between states - no peptide correspondence
#' required at all. Useful when peptide identification differs enough
#' between states (chimeras, point mutants, differently-digested runs) that
#' peptide-matched comparison would discard too much data, at the cost of a
#' real confound: two states can appear to differ purely because their
#' peptide *maps* differ (different lengths/boundaries), not because the
#' underlying exchange does. See Details.
#'
#' @section This is exploratory, not confirmatory:
#' Unlike [compare_states_tp()], results here are `shift_hl`/`ci_low`/
#' `ci_high`/`p_value` - deliberately not `is_sig`/`delta` - so this
#' function's output is never mistaken for the package's confirmatory
#' Hageman & Weis criterion. Always inspect the returned `coverage` table
#' before trusting a result: if `frac_shared_all_states` is low for a
#' region, the states' peptide maps differ enough there that some of any
#' observed shift may reflect peptide composition, not biology.
#'
#' @section Why `regions` is required:
#' Pooling every peptide in the whole protein mixes a strongly bimodal mix
#' of protected-core and exposed-loop behavior; the pooled median is then
#' governed mostly by which peptides happened to be identified (composition
#' has more sequence available in exposed regions to fully deuterate),
#' not by biology. Restricting to one structurally coherent region at a
#' time keeps the pooled distribution closer to unimodal and interpretable.
#'
#' @section Read the per-peptide diagnostic first:
#' `summary` lists the matched per-peptide diagnostic first
#' (`n_peptides_tested`, `n_peptides_significant`,
#' `frac_peptides_significant`), then the pooled region test. The pooled
#' test is approximate: it treats the region's peptides as independent
#' samples although overlapping peptides share residues, and it compares
#' peptide means without their replicate variance. For peptides measured in
#' every state, [compare_regions()] (a two-way ANOVA with peptide as a
#' block) is the more powerful region-level test.
#'
#' @section Pooled significance is not peptide significance:
#' The pooled test above answers "does this region's *typical* exchange
#' shift" - it says nothing about whether *some* peptides in the region
#' shift while others don't. A real, spatially localized effect (e.g. a
#' binding interface covering a handful of peptides inside an otherwise
#' unaffected 200-residue domain) gets averaged against every unaffected
#' peptide riding along in the same region and can come back completely
#' non-significant even though individual peptides in it clearly differ.
#' To surface this, `summary` also reports `n_peptides_tested`,
#' `n_peptides_significant`, and `frac_peptides_significant`: for every
#' peptide identity (Start, End, Charge and Sequence - each charge state
#' counts as its own peptide) shared between that pair's two states in
#' that region, a
#' per-peptide matched Welch t-test (replicates of that one peptide, state1
#' vs. state2 - the same kind of test [compare_states_tp()] and
#' `plot_fissure_vent()` use), with p-values adjusted together across every
#' peptide-by-pair test in this call (same family as `p_adj` above) and
#' called significant at `p_adj < pvalue`. A region with a flat/non-significant
#' pooled result but `frac_peptides_significant > 0` means real, localized
#' signal exists in there - shrink `regions` to isolate it, or use
#' peptide-matched testing directly instead of pooling.
#'
#' @section What gets compared: `value`:
#' `filepath`'s raw `# Deut` column (Da uptake, read as `X..Deut`) is used
#' directly by default (`value = "raw"`) - no per-peptide size correction of
#' any kind. Two other options exist, off by default:
#' \describe{
#'   \item{`"raw"`}{Default. Raw deuterium uptake in Da, exactly the `# Deut`
#'     column from the input file - the most direct reading of the data, and
#'     what `compare_states_tp()`/`plot_fissure_vent()`/`diff_peptides_prep()`
#'     use as well (their `percent = FALSE` default).}
#'   \item{`"percent"`}{HDX Examiner's own `Deut %` column - uptake already
#'     expressed relative to that peptide's theoretical maximum.}
#'   \item{`"frac_deut"`}{This package's own derived metric: `# Deut`
#'     divided by the peptide's number of exchangeable protons (see
#'     [nb_exch_deut()]) - a finer-grained per-residue normalization than
#'     `"percent"`. Useful specifically when pooling peptides of very
#'     different lengths together in one region and the *length* difference
#'     itself would otherwise dominate the pooled comparison (a longer
#'     peptide mechanically carries a larger raw Da value regardless of its
#'     intrinsic exchange rate) - trades that off against normalizing away
#'     part of the signal, which is why it is opt-in, not the default.}
#' }
#' The per-peptide p-values do not depend on this choice (dividing one
#' peptide's two states by the same constant does not change a t-test's
#' p-value), but the critical-interval half of the per-peptide rule and the
#' *pooled* result (`shift_hl`/`p_value`/`p_omnibus`, which pools values
#' across different peptides) are computed on the chosen scale.
#'
#' @section `method`:
#' \describe{
#'   \item{`"peptide"` (default)}{Each peptide contributes one value per
#'     state (its replicates averaged first). Recommended default -
#'     avoids replicate pseudo-replication (replicates of one peptide are
#'     not independent draws of "how this state exchanges"). Omnibus:
#'     Kruskal-Wallis. Pairwise: Wilcoxon rank-sum, reported alongside the
#'     Hodges-Lehmann shift estimate.}
#'   \item{`"classic"`}{Pools every replicate row as its own observation
#'     (`replicates` x n_peptides points per state) - reproduces the
#'     pooling this package's own prior ad hoc analyses used. Triggers a
#'     warning: replicates within a peptide are highly correlated, so this
#'     inflates the effective sample size and understates p-values. Kept
#'     for reproducing prior results, not recommended for new analyses.}
#'   \item{`"resample"`}{Peptide-level like `"peptide"`, but p-values come
#'     from a label-permutation test (permuting which peptide belongs to
#'     which state and recomputing the test statistic) rather than the
#'     Kruskal-Wallis/Wilcoxon asymptotic approximation, and the pairwise
#'     `shift_hl` gets a percentile bootstrap confidence interval
#'     (resampling peptides with replacement within each state). Slower
#'     (`n_perm` + `n_boot` resamples per pair/region); more defensible
#'     with small peptide counts per state where the asymptotic tests are
#'     unreliable.}
#' }
#'
#' @param filepath filepath to the input file (read via
#'   [all_data_long_input()] - no peptide matching between states happens
#'   at any stage).
#' @param times a single `Deut.Time` value to compare at.
#' @param regions required named list of `c(start, end)` residue ranges,
#'   e.g. `list(NTR = c(1, 85), ACD = c(86, 165))` - see Details for why
#'   whole-protein pooling isn't offered.
#' @param states character vector of states to compare. Default `NULL`
#'   uses every state found in `filepath`.
#' @param value one of `"raw"` (default), `"percent"`, `"frac_deut"` - which
#'   column drives the comparison. See "What gets compared: `value`".
#' @param method one of `"peptide"` (default), `"classic"`, `"resample"` -
#'   see Details.
#' @param reference optional state name - if given, only `reference` vs.
#'   each other state is tested, instead of every pair.
#' @param replicates number of replicates per state. Default `3`.
#' @param pvalue significance cutoff. This function does not compute
#'   `is_sig` for the pooled `shift_hl`/`p_value` result itself (callers can
#'   filter `p_adj < pvalue` there themselves), but `pvalue` does gate the
#'   per-peptide `frac_peptides_significant` diagnostic - see "Pooled
#'   significance is not peptide significance". Default `0.01`.
#' @param alpha alpha of the pooled Hageman & Weis critical interval used by
#'   the per-peptide diagnostic. Default `0.01`.
#' @param p_adjust multiple-testing correction passed to [stats::p.adjust()]:
#'   `"none"` (default) or e.g. `"BH"`, applied across every region-by-pair
#'   test together (and separately across every region-by-pair-by-peptide
#'   test for the per-peptide diagnostic).
#' @param p_adjust_method deprecated; earlier versions used
#'   `p_adjust_method = "BH"` by default. Supplying it sets `p_adjust`, with
#'   a warning.
#' @param seed optional random seed for `method = "resample"`, making the
#'   permutation p-values and bootstrap interval reproducible. The caller's
#'   own random-number stream is left unchanged. Default `NULL`.
#' @param n_perm number of label permutations for `method = "resample"`.
#'   Default `1000`.
#' @param n_boot number of bootstrap resamples for `method = "resample"`'s
#'   confidence interval. Default `1000`.
#' @return a list with `value` (the `value` argument actually used, echoed
#'   back for plotting/labeling) and three data.frames: `summary` (one row
#'   per region-by-state-pair: `region`, `state1`, `state2`, `n1`, `n2`,
#'   `shift_hl`, `ci_low`, `ci_high` - `NA` unless `method = "resample"` -
#'   `p_value`, `p_omnibus`, `p_adj`, `n_peptides_tested`,
#'   `n_peptides_significant`, `frac_peptides_significant` - see "Pooled
#'   significance is not peptide significance"), `coverage` (one row per
#'   region-by-state: `n_peptides`, `mean_length`, and
#'   `frac_shared_all_states` - the fraction of that region's peptide
#'   identities, by Start/End/Sequence, present in *every* state being
#'   compared, a diagnostic for how much the matched-vs-unmatched confound
#'   applies here), and `data` (the pooled `deut_value` values the tests
#'   actually ran on - one row per region-by-state-by-peptide, or per
#'   replicate under `method = "classic"` - for plotting with
#'   [plot_compare_states_unmatched()]), and `peptides` (the matched
#'   per-peptide results behind the diagnostic: `delta`, `p_value`, `p_adj`,
#'   `ci` and `significant` per region, pair and peptide).
#' @seealso [compare_states_tp()], [all_data_long_input()],
#'   [residue_level_average()], [plot_compare_states_unmatched()],
#'   [plot_coverage_diagnostics()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- compare_states_unmatched(file_nm, times = "1800.00s",
#'                                  regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' res$summary
#' res$coverage
#' }
#' @export
compare_states_unmatched <- function(filepath, times, regions, states = NULL,
                                      value = c("raw", "percent", "frac_deut"),
                                      method = c("peptide", "classic", "resample"),
                                      reference = NULL,
                                      replicates = 3, pvalue = 0.01, alpha = 0.01,
                                      p_adjust = "none", n_perm = 1000, n_boot = 1000, seed = NULL,
                                      p_adjust_method) {
  value <- match.arg(value)
  method <- match.arg(method)
  if (!missing(p_adjust_method)) {
    warning("p_adjust_method is deprecated: use p_adjust (default now \"none\").", call. = FALSE)
    p_adjust <- p_adjust_method
  }
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  if (!is.null(seed)) {
    ## reproducible resampling without changing the caller's random stream
    old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }
  if (missing(regions) || is.null(regions) || length(regions) == 0) {
    stop("regions is required - whole-protein pooling mixes very different exchange behavior ",
         "(protected core vs. exposed loops) and is not a reliable summary. Supply at least one ",
         "named c(start, end) range, e.g. regions = list(my_domain = c(10, 90)).")
  }

  b <- all_data_long_input(filepath)
  b <- b[b$Deut.Time == times, ]
  if (nrow(b) == 0) {
    stop(sprintf("no rows with Deut.Time == '%s'", times))
  }
  b$deut_value <- switch(value,
    raw = as.numeric(b$X..Deut),
    percent = as.numeric(b$Deut..),
    frac_deut = b$frac_deut
  )
  if (is.null(states)) states <- unique(b$Protein.State)
  if (length(states) < 2) {
    stop("need at least two states to compare.")
  }
  if (!is.null(reference) && !(reference %in% states)) {
    stop(sprintf("reference '%s' not among states: %s", reference, paste(states, collapse = ", ")))
  }
  b <- b[b$Protein.State %in% states, ]
  ## Charge is part of the identity: each charge state is its own peptide
  ## (as in output_tp()/compare_states_tp()), so a peptide seen at two charge
  ## states is never pooled as 2x the replicates in the matched t-test
  b$peptide_id <- paste(b$Start, b$End, paste0("z", b$Charge), b$Sequence, sep = "_")

  if (method == "classic") {
    warning("method = \"classic\" pools every replicate row as its own observation - replicates ",
             "within a peptide are highly correlated (not independent draws), which inflates the ",
             "effective sample size and understates p-values. method = \"peptide\" (the default) ",
             "avoids this; use \"classic\" only to reproduce a prior analysis built on raw-replicate pooling.",
             call. = FALSE)
  }

  pairs <- if (!is.null(reference)) {
    lapply(setdiff(states, reference), function(s) c(reference, s))
  } else {
    utils::combn(states, 2, simplify = FALSE)
  }

  ## one value per peptide per state (replicates averaged), unless
  ## "classic" pools raw replicate rows directly
  if (method == "classic") {
    unit_df <- b[, c("Protein.State", "Start", "End", "peptide_id", "deut_value")]
  } else {
    unit_df <- stats::aggregate(deut_value ~ Protein.State + Start + End + peptide_id, data = b, FUN = mean)
  }

  results <- list(); coverage <- list(); pooled_data <- list(); peptide_sig_list <- list()
  for (rname in names(regions)) {
    region <- regions[[rname]]
    in_r <- unit_df$Start >= region[1] & unit_df$End <= region[2]
    ud <- unit_df[in_r, ]
    if (nrow(ud) > 0) {
      pooled_data[[rname]] <- data.frame(region = rname, ud, stringsAsFactors = FALSE)
    }
    b_region <- b[b$Start >= region[1] & b$End <= region[2], ]

    ids_by_state <- lapply(states, function(s) unique(ud$peptide_id[ud$Protein.State == s]))
    shared <- Reduce(intersect, ids_by_state)
    all_ids <- Reduce(union, ids_by_state)
    frac_shared <- if (length(all_ids) > 0) length(shared) / length(all_ids) else NA_real_

    coverage[[rname]] <- do.call(rbind, lapply(states, function(s) {
      sd_ <- ud[ud$Protein.State == s, ]
      data.frame(region = rname, state = s,
                 n_peptides = length(unique(sd_$peptide_id)),
                 mean_length = mean(sd_$End - sd_$Start + 1, na.rm = TRUE),
                 frac_shared_all_states = frac_shared,
                 stringsAsFactors = FALSE)
    }))

    if (nrow(ud) == 0) next

    p_omnibus <- NA_real_
    if (length(unique(ud$Protein.State)) >= 2 && length(unique(ud$deut_value)) >= 2) {
      if (method == "resample") {
        obs_stat <- tryCatch(stats::kruskal.test(ud$deut_value, factor(ud$Protein.State))$statistic,
                              error = function(e) NA_real_)
        if (!is.na(obs_stat)) {
          perm_stat <- vapply(seq_len(n_perm), function(i) {
            tryCatch(stats::kruskal.test(ud$deut_value, factor(sample(ud$Protein.State)))$statistic,
                      error = function(e) NA_real_)
          }, numeric(1))
          p_omnibus <- (1 + sum(perm_stat >= obs_stat, na.rm = TRUE)) / (1 + n_perm)
        }
      } else {
        p_omnibus <- tryCatch(stats::kruskal.test(ud$deut_value, factor(ud$Protein.State))$p.value,
                               error = function(e) NA_real_)
      }
    }

    pair_rows <- lapply(pairs, function(pr) {
      s1 <- pr[1]; s2 <- pr[2]
      x <- ud$deut_value[ud$Protein.State == s1]
      y <- ud$deut_value[ud$Protein.State == s2]
      if (length(x) < 2 || length(y) < 2) {
        return(data.frame(region = rname, state1 = s1, state2 = s2, n1 = length(x), n2 = length(y),
                           shift_hl = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                           p_value = NA_real_, p_omnibus = p_omnibus, stringsAsFactors = FALSE))
      }
      shift <- .hl_shift(x, y)
      ci_low <- NA_real_; ci_high <- NA_real_
      if (method == "resample") {
        p_val <- .perm_test_hl(x, y, n_perm)
        boot_shifts <- vapply(seq_len(n_boot), function(i) {
          .hl_shift(sample(x, length(x), replace = TRUE), sample(y, length(y), replace = TRUE))
        }, numeric(1))
        qs <- stats::quantile(boot_shifts, c(0.025, 0.975), na.rm = TRUE)
        ci_low <- qs[[1]]; ci_high <- qs[[2]]
      } else {
        p_val <- tryCatch(stats::wilcox.test(y, x, exact = FALSE)$p.value, error = function(e) NA_real_)
      }
      pep_pv <- .peptide_matched_pvals(b_region, s1, s2)
      if (nrow(pep_pv) > 0) {
        ## pooled Hageman & Weis CI from the matched peptides' replicate SDs
        pep_pv$ci <- CI_2pts(pep_pv$sd1, pep_pv$sd2, alpha = alpha, n1 = pep_pv$n1, n2 = pep_pv$n2)
        peptide_sig_list[[length(peptide_sig_list) + 1]] <<- data.frame(
          region = rname, state1 = s1, state2 = s2, pep_pv, stringsAsFactors = FALSE)
      }

      data.frame(region = rname, state1 = s1, state2 = s2, n1 = length(x), n2 = length(y),
                 shift_hl = shift, ci_low = ci_low, ci_high = ci_high,
                 p_value = p_val, p_omnibus = p_omnibus, stringsAsFactors = FALSE)
    })
    results[[rname]] <- do.call(rbind, pair_rows)
  }

  summary_df <- do.call(rbind, results)
  if (is.null(summary_df) || nrow(summary_df) == 0) {
    stop("no peptides found in any requested region.")
  }
  summary_df$p_adj <- stats::p.adjust(summary_df$p_value, method = p_adjust)

  ## per-peptide matched significance rate per region-by-pair: catches a
  ## real but spatially localized effect that the pooled test above dilutes
  ## away (see @section "Pooled significance is not peptide significance").
  ## Same rule as the rest of the package: Welch p (raw, or adjusted with
  ## p_adjust across every peptide-by-pair test in this call) < pvalue AND
  ## |delta| > the pooled critical interval of that region and pair.
  peptide_sig_df <- do.call(rbind, peptide_sig_list)
  summary_df$n_peptides_tested <- 0L
  summary_df$n_peptides_significant <- 0L
  summary_df$frac_peptides_significant <- NA_real_
  if (!is.null(peptide_sig_df) && nrow(peptide_sig_df) > 0) {
    peptide_sig_df$p_adj <- stats::p.adjust(peptide_sig_df$p_value, method = p_adjust)
    peptide_sig_df$significant <- !is.na(peptide_sig_df$p_adj) & peptide_sig_df$p_adj < pvalue &
      !is.na(peptide_sig_df$delta) & abs(peptide_sig_df$delta) > peptide_sig_df$ci
    key_pep <- paste(peptide_sig_df$region, peptide_sig_df$state1, peptide_sig_df$state2, sep = "||")
    key_sum <- paste(summary_df$region, summary_df$state1, summary_df$state2, sep = "||")
    n_tested <- tapply(!is.na(peptide_sig_df$p_value), key_pep, sum)
    n_sig <- tapply(peptide_sig_df$significant, key_pep, sum, na.rm = TRUE)
    m <- match(key_sum, names(n_tested))
    summary_df$n_peptides_tested <- ifelse(is.na(m), 0L, n_tested[m])
    summary_df$n_peptides_significant <- ifelse(is.na(m), 0L, n_sig[m])
    summary_df$frac_peptides_significant <- ifelse(summary_df$n_peptides_tested > 0,
                                                     summary_df$n_peptides_significant / summary_df$n_peptides_tested,
                                                     NA_real_)
  }
  ## the matched per-peptide diagnostic first; the pooled (approximate)
  ## region test after it
  first <- c("region", "state1", "state2", "n_peptides_tested", "n_peptides_significant",
             "frac_peptides_significant")
  summary_df <- summary_df[, c(first, setdiff(colnames(summary_df), first))]
  rownames(summary_df) <- NULL
  if (!is.null(peptide_sig_df)) rownames(peptide_sig_df) <- NULL

  coverage_df <- do.call(rbind, coverage)
  rownames(coverage_df) <- NULL

  data_df <- do.call(rbind, pooled_data)
  rownames(data_df) <- NULL

  .hdx_stamp(list(summary = summary_df, coverage = coverage_df, data = data_df, peptides = peptide_sig_df,
                  value = value),
             list(times = times, regions = regions, states = states, value = value, method = method,
                  pvalue = pvalue, alpha = alpha, p_adjust = p_adjust, n_perm = n_perm, n_boot = n_boot,
                  seed = if (is.null(seed)) "none" else seed))
}
