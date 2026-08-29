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
#' @param method one of `"peptide"` (default), `"classic"`, `"resample"` -
#'   see Details.
#' @param reference optional state name - if given, only `reference` vs.
#'   each other state is tested, instead of every pair.
#' @param replicates number of replicates per state. Default `3`.
#' @param pvalue significance cutoff - informational only here (this
#'   function does not compute `is_sig`); provided so callers can filter
#'   `p_adj < pvalue` themselves. Default `0.01`.
#' @param p_adjust_method passed to [stats::p.adjust()], applied across
#'   every region-by-pair test together. Default `"BH"`.
#' @param n_perm number of label permutations for `method = "resample"`.
#'   Default `1000`.
#' @param n_boot number of bootstrap resamples for `method = "resample"`'s
#'   confidence interval. Default `1000`.
#' @return a list with three data.frames: `summary` (one row per
#'   region-by-state-pair: `region`, `state1`, `state2`, `n1`, `n2`,
#'   `shift_hl`, `ci_low`, `ci_high` - `NA` unless `method = "resample"` -
#'   `p_value`, `p_omnibus`, `p_adj`), `coverage` (one row per
#'   region-by-state: `n_peptides`, `mean_length`, and
#'   `frac_shared_all_states` - the fraction of that region's peptide
#'   identities, by Start/End/Sequence, present in *every* state being
#'   compared, a diagnostic for how much the matched-vs-unmatched confound
#'   applies here), and `data` (the pooled values the tests actually ran
#'   on - one row per region-by-state-by-peptide, or per replicate under
#'   `method = "classic"` - for plotting with [plot_compare_states_unmatched()]).
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
                                      method = c("peptide", "classic", "resample"),
                                      reference = NULL,
                                      replicates = 3, pvalue = 0.01,
                                      p_adjust_method = "BH", n_perm = 1000, n_boot = 1000) {
  method <- match.arg(method)
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
  if (is.null(states)) states <- unique(b$Protein.State)
  if (length(states) < 2) {
    stop("need at least two states to compare.")
  }
  if (!is.null(reference) && !(reference %in% states)) {
    stop(sprintf("reference '%s' not among states: %s", reference, paste(states, collapse = ", ")))
  }
  b <- b[b$Protein.State %in% states, ]
  b$peptide_id <- paste(b$Start, b$End, b$Sequence, sep = "_")

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
    unit_df <- b[, c("Protein.State", "Start", "End", "peptide_id", "frac_deut")]
  } else {
    unit_df <- stats::aggregate(frac_deut ~ Protein.State + Start + End + peptide_id, data = b, FUN = mean)
  }

  results <- list(); coverage <- list(); pooled_data <- list()
  for (rname in names(regions)) {
    region <- regions[[rname]]
    in_r <- unit_df$Start >= region[1] & unit_df$End <= region[2]
    ud <- unit_df[in_r, ]
    if (nrow(ud) > 0) {
      pooled_data[[rname]] <- data.frame(region = rname, ud, stringsAsFactors = FALSE)
    }

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
    if (length(unique(ud$Protein.State)) >= 2 && length(unique(ud$frac_deut)) >= 2) {
      if (method == "resample") {
        obs_stat <- tryCatch(stats::kruskal.test(ud$frac_deut, factor(ud$Protein.State))$statistic,
                              error = function(e) NA_real_)
        if (!is.na(obs_stat)) {
          perm_stat <- vapply(seq_len(n_perm), function(i) {
            tryCatch(stats::kruskal.test(ud$frac_deut, factor(sample(ud$Protein.State)))$statistic,
                      error = function(e) NA_real_)
          }, numeric(1))
          p_omnibus <- (1 + sum(perm_stat >= obs_stat, na.rm = TRUE)) / (1 + n_perm)
        }
      } else {
        p_omnibus <- tryCatch(stats::kruskal.test(ud$frac_deut, factor(ud$Protein.State))$p.value,
                               error = function(e) NA_real_)
      }
    }

    pair_rows <- lapply(pairs, function(pr) {
      s1 <- pr[1]; s2 <- pr[2]
      x <- ud$frac_deut[ud$Protein.State == s1]
      y <- ud$frac_deut[ud$Protein.State == s2]
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
  summary_df$p_adj <- stats::p.adjust(summary_df$p_value, method = p_adjust_method)
  rownames(summary_df) <- NULL

  coverage_df <- do.call(rbind, coverage)
  rownames(coverage_df) <- NULL

  data_df <- do.call(rbind, pooled_data)
  rownames(data_df) <- NULL

  list(summary = summary_df, coverage = coverage_df, data = data_df)
}
