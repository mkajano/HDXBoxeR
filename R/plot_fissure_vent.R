#' Fissure_vent plot
#'
#' modification of the volcano plot.
#' The volcano plot linearized.
#' x and y axis are the deuteration for two states
#' Deuteration can be as either % deuteration, uptake, or fraction deuteration
#'
#' The CI is plotted as the lines parallel to the diagonal
#' The uncerainity of the peptide deuteration is plotted as a ellipsoid using the pvalues.
#' The pvalues were calculated using the Welch ttest and they are converted to the SD using the equation:
#' SD=sqrt(replicates)*(uper_limit-lower_limit)/t
#' https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
#'
#' The function as a default compares the defined states.
#' To be able to compare more states use loop.
#'
#' Peptides are drawn with a solid black outline when they clear BOTH the
#' hybrid significance criteria (outside the global CI band AND Welch
#' p < pvalue) - the same two-part standard as the rest of this package
#' (CI_2pts + a p-value cutoff). Peptides clearing only one, or neither,
#' are drawn in grey - visible, but not flagged as significant.
#'
#' @param df output of output_tp() function.
#' @param times Deut times
#' @param pvalue cutoff for the pvalue. A peptide is only outlined as
#'   significant if BOTH abs(state2-state1) exceeds the CI band AND its
#'   Welch t-test p-value (state1 vs state2, via pv_timepoint) is below
#'   this cutoff.
#' @param replicates number of replicates in sample. Default set to 3.
#' @param state1 First protein state. Matched by substring against the
#'   averaged/SD column names (e.g. "Fresh" matches "av_Fresh_pH7"), so it
#'   does not need to be the exact column name - just specific enough to
#'   match exactly one column. Defaults to the first state found in `df`.
#' @param state2 Second protein state, same matching rule. Defaults to the
#'   second state found in `df`.
#' @param sd_source flag. "CI" (default) calculates SD from the critical
#'   interval (sd_ci_timepoint); anything else uses the raw SD of the
#'   peptide's replicates (sd_timepoint).
#' @param ylim y axis limit
#' @param xlim x axis limit
#' @param ... other plotting parameters, passed to DrawEllipse
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @return Fissure vent plot showing peptide-level deuteration in state1 vs
#'   state2, with the CI band and per-peptide uncertainty ellipses.
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' plot_fissure_vent(df = a, replicates = 3, times = 3.00, alpha = 0.01,
#'                    pvalue = 0.01, state1 = "Unbound", state2 = "bound")
#' }
#' @export

plot_fissure_vent <- function(df, replicates = 3, times, state1, state2, alpha = 0.01, pvalue = 0.01,
                               position_col = NA, deut = "", sd_source = "CI", length1 = NA,
                               xlim = NULL, ylim = NULL, ...) {

  av1 <- ave_timepoint(df, replicates)

  nm1 <- str_sub(colnames(av1)[7:ncol(av1)], start = 4, end = -10)
  if (missing(state1)) state1 <- nm1[1]
  if (missing(state2)) state2 <- nm1[2]

  # av1's columns are prefixed "av_...", sd1's are "sd_..."/"sd_ci_..." - the
  # SAME state string can't index both directly (this was the original bug:
  # state1/state2 only ever worked as a raw column name for ONE of the two
  # data frames). Resolve by matching state against the START of the root
  # name (after stripping any known prefix), not just "contains" - a plain
  # "contains" match is ambiguous whenever one state name is a substring of
  # another, e.g. this package's own shipped example has states "Unbound"
  # and "bound" ("bound" is a substring of "Unbound") - caught by testing
  # this fix against that real example data, not assumed safe.
  resolve_col <- function(state, cols, known_prefixes = character(0)) {
    root <- cols
    for (p in known_prefixes) root <- sub(paste0("^", p), "", root)
    hits <- cols[startsWith(root, state)]
    if (length(hits) == 0) {
      stop(sprintf("No column matching state '%s' found among: %s", state, paste(cols, collapse = ", ")))
    }
    if (length(hits) > 1) {
      stop(sprintf("State '%s' matches multiple columns (%s) - use a more specific name",
                    state, paste(hits, collapse = ", ")))
    }
    hits
  }

  if (sd_source == "CI") {
    sd1 <- sd_ci_timepoint(df, replicates, alpha)
  } else {
    sd1 <- sd_timepoint(df, replicates)
  }

  av_state1 <- resolve_col(state1, colnames(av1), "av_")
  av_state2 <- resolve_col(state2, colnames(av1), "av_")
  sd_prefixes <- c("sd_ci_", "sd_")  # check the longer prefix first
  sd_state1 <- resolve_col(state1, colnames(sd1), sd_prefixes)
  sd_state2 <- resolve_col(state2, colnames(sd1), sd_prefixes)

  # Per-peptide Welch p-value for state1 vs state2 specifically, computed
  # directly from df's raw replicate columns - NOT via pv_timepoint(), which
  # always compares every state against whichever state sits in df's column
  # 7, regardless of what state1/state2 are here. That mismatch would
  # silently mislabel significance whenever state1 isn't df's first state.
  resolve_raw_cols <- function(state, cols, replicates) {
    hits <- cols[startsWith(cols, state)]
    if (length(hits) != replicates) {
      stop(sprintf("Expected %d raw replicate columns matching state '%s', found %d: %s",
                    replicates, state, length(hits), paste(hits, collapse = ", ")))
    }
    hits
  }
  raw_cols1 <- resolve_raw_cols(state1, colnames(df), replicates)
  raw_cols2 <- resolve_raw_cols(state2, colnames(df), replicates)
  pv_state1v2 <- vapply(seq_len(nrow(df)), function(r) {
    tryCatch(t.test(as.numeric(df[r, raw_cols1]), as.numeric(df[r, raw_cols2]))$p.value,
             error = function(e) NA_real_)
  }, numeric(1))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(oma = c(3, 3, 1, 1))

  plot(av1[, av_state1], av1[, av_state1], type = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim, ...)
  mtext("Deuteration (state 1)", side = 1, line = 0.7, outer = TRUE)
  mtext("Deuteration (state 2)", side = 2, line = 0.7, outer = TRUE)
  points(c(-10, 110), c(-10, 110), type = "l")

  draw_one <- function(idx, ci, col_default, use_position_col) {
    n_skipped <- 0
    for (pep in idx) {
      # a peptide can have an NA SD (e.g. zero-variance replicates make the
      # underlying t-test undefined) - skip drawing its ellipse rather than
      # crashing DrawEllipse/the significance check on a missing value,
      # caught by testing this against real replicate data, not assumed away
      if (anyNA(c(av1[pep, av_state1], av1[pep, av_state2], sd1[pep, sd_state1], sd1[pep, sd_state2]))) {
        n_skipped <- n_skipped + 1
        next
      }
      delta <- abs(av1[pep, av_state2] - av1[pep, av_state1])
      is_sig <- (delta > ci / 2) && !is.na(pv_state1v2[pep]) && (pv_state1v2[pep] < pvalue)
      border_col <- if (use_position_col) position_col[pep] else if (is_sig) "black" else "grey60"
      fill_col   <- if (use_position_col) position_col[pep] else if (is_sig) "grey30" else "grey85"
      DrawEllipse(av1[pep, av_state1], av1[pep, av_state2],
                  sd1[pep, sd_state1] / 2, sd1[pep, sd_state2] / 2,
                  col = fill_col, border = border_col, ...)
    }
    if (n_skipped > 0) message(sprintf("  (%d peptide(s) skipped: NA average or SD, e.g. zero-variance replicates)", n_skipped))
  }

  use_position_col <- !is.na(position_col[1])
  for (i in seq_along(times)) {
    ind1 <- select_indices(av1, times = times[i], length = length1)
    ci <- CI_2pts(stats::na.omit(sd1[ind1, sd_state1]), stats::na.omit(sd1[ind1, sd_state2]), replicates, alpha)
    message(sprintf("time %s: CI = %.3f", times[i], ci))
    points(c(-10, 110), c(-10, 110) + ci / 2, type = "l", col = "grey")
    points(c(-10, 110), c(-10, 110) - ci / 2, type = "l", col = "grey")
    draw_one(ind1, ci, "grey30", use_position_col)
  }
}
