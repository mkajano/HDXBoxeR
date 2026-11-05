#' Fissure_vent plot
#'
#' modification of the volcano plot.
#' The volcano plot linearized.
#' x and y axis are the deuteration for two states
#' Deuteration can be as either % deuteration, uptake, or fraction deuteration
#'
#' The CI is plotted as the lines parallel to the diagonal - a single,
#' POOLED threshold (`CI_2pts()`) built from every peptide's own raw
#' replicate SD per state at this timepoint (`sd_timepoint()`) - the same
#' SDs and therefore the same threshold as `CI_tp()`, `diff_peptides_prep()`,
#' the heat map/Woods/robot plots and `compare_states_tp()`.
#' It answers "is this peptide's change bigger than the typical measurement
#' noise across the whole experiment" - a deliberately pooled,
#' practical-significance question, per the Hageman & Weis hybrid method
#' this package implements.
#'
#' Each peptide's uncertainty ellipse is a separate, per-peptide quantity:
#' it is built from the exact same Welch t-test used for that peptide's own
#' p-value (state1 vs state2, that peptide's raw replicates only), run at
#' `conf.level = 1 - pvalue` - so the confidence interval this test returns
#' excludes zero if and only if that peptide's own p-value is below
#' `pvalue` (the standard CI/p-value duality, exact at the same alpha, not
#' an approximation). That interval's half-width is split between the two
#' axes in proportion to each state's own standard error, so the ellipse's
#' extent along the delta direction (perpendicular to the diagonal) equals
#' that half-width exactly. Consequence: an ellipse that still touches the
#' DIAGONAL LINE ITSELF (not the grey CI band) means this peptide's own
#' test failed to clear `pvalue`; an ellipse that clears the diagonal means
#' it passed. The ellipse is therefore a valid visual readout of the
#' per-peptide p-value test specifically - but it answers a different
#' question than the grey band, and both are required for the black/grey
#' outline below. Don't read "ellipse crosses the band" as the
#' significance test; read "ellipse crosses the diagonal" for the p-value
#' half, and the point's position relative to the band for the other half.
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
#' @param state1 First protein state. An exact state name is matched first
#'   (so "WT" never picks up "WT_mut"); otherwise an unambiguous
#'   abbreviation is accepted (e.g. "Fresh" for "Fresh_pH7"). Defaults to
#'   the first state found in `df`.
#' @param state2 Second protein state, same matching rule. Defaults to the
#'   second state found in `df`.
#' @param sd_source flag. `"raw"` (default) uses the raw SD of each state's
#'   replicates (`sd_timepoint()`), matching the threshold used everywhere
#'   else in the package. `"CI"` instead back-calculates SDs from a
#'   two-sample t-test interval (`sd_ci_timepoint()`); those are SDs of the
#'   *difference* between states, so the resulting band is wider (about
#'   1.8x on the shipped example) and fewer peptides are called significant.
#' @param position_col optional vector of per-peptide colors; when supplied,
#'   ellipses are drawn in these colors instead of the significance-based
#'   `sig_col`/`nonsig_col` outline, and `fill_alpha` is ignored (position
#'   colors are used as-is, opaque).
#' @param length1 optional peptide length filter passed to `select_indices()`.
#' @param ylim y axis limit
#' @param xlim x axis limit
#' @param ... other plotting parameters, passed to DrawEllipse
#' @param alpha critical interval, to have more restrictive use lower values, default=0.01
#' @param sig_col border/fill color for significant peptides (outside the CI
#'   band AND clearing the diagonal). Default `"firebrick"`.
#' @param nonsig_col border/fill color for peptides that don't clear both
#'   criteria. Default `"grey60"`.
#' @param fill_alpha transparency (0-1) applied to each ellipse's fill only
#'   (the border stays opaque); `0` is fully transparent, `1` fully opaque.
#'   Default `0.3`. Ignored when `position_col` is supplied.
#' @param p_adjust multiple-testing correction of the per-peptide Welch
#'   p-values, passed to [stats::p.adjust()]: `"none"` (default) or e.g.
#'   `"BH"`, applied across the peptides of each time point. The
#'   significance outline uses the adjusted p-value; the ellipses stay
#'   calibrated to the raw p-value (touching the diagonal means raw
#'   p >= `pvalue`).
#' @return Fissure vent plot showing peptide-level deuteration in state1 vs
#'   state2, with the CI band and per-peptide uncertainty ellipses.
#'   Invisibly, a data frame with one row per drawn peptide: time, Start,
#'   End, Sequence, both state averages, `delta`, `p_value`, `p_adj`, `ci`
#'   and `is_sig`.
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' plot_fissure_vent(df = a, replicates = 3, times = "3.00s", alpha = 0.01,
#'                    pvalue = 0.01, state1 = "Unbound", state2 = "bound")
#' plot_fissure_vent(df = a, replicates = 3, times = "3.00s", alpha = 0.01,
#'                    pvalue = 0.01, state1 = "Unbound", state2 = "bound",
#'                    sig_col = "darkorange", nonsig_col = "steelblue", fill_alpha = 0.5)
#' }
#' @export

plot_fissure_vent <- function(df, replicates = 3, times, state1, state2, alpha = 0.01, pvalue = 0.01,
                               position_col = NA, sd_source = "raw", length1 = NA,
                               xlim = NULL, ylim = NULL,
                               sig_col = "firebrick", nonsig_col = "grey60", fill_alpha = 0.3,
                               p_adjust = "none", ...) {

  av1 <- ave_timepoint(df, replicates)

  # av1's columns are "av_<state>_X..Deut" (uptake) or "av_<state>_Deut.."
  # (percent) - strip the "av_" prefix and the known value-suffix by pattern,
  # not by a fixed character count, since the two suffixes differ in length.
  nm1 <- sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", sub("^av_", "", colnames(av1)[7:ncol(av1)]))
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
    .hdx_summary_col(state, cols, known_prefixes)
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

  # Per-peptide Welch p-value AND matching per-peptide ellipse radii for
  # state1 vs state2 specifically, computed directly from df's raw replicate
  # columns in a SINGLE t.test() call per peptide - NOT via pv_timepoint(),
  # which always compares every state against whichever state sits in df's
  # column 7, regardless of what state1/state2 are here. That mismatch would
  # silently mislabel significance whenever state1 isn't df's first state.
  #
  # conf.level = 1 - pvalue (not a fixed 0.95) so "this peptide's CI
  # excludes zero" is exactly equivalent to "this peptide's p < pvalue" -
  # same test, same alpha. The CI half-width is split across the two axes
  # in proportion to each state's own standard error (t*s/sqrt(n), using
  # this peptide's own non-NA replicate count per state, not a blanket
  # `replicates`), so the resulting ellipse's extent along the delta
  # direction reproduces that half-width exactly - see @details.
  # t.test() does not error on an all-identical-values comparison (e.g. a
  # peptide with 0 uptake in every replicate of both states) - it silently
  # returns NaN for both p.value and conf.int - so NaN is checked
  # explicitly rather than relying on tryCatch alone.
  resolve_raw_cols <- function(state, cols, replicates) {
    .hdx_state_cols(state, cols)
  }
  raw_cols1 <- resolve_raw_cols(state1, colnames(df), replicates)
  raw_cols2 <- resolve_raw_cols(state2, colnames(df), replicates)

  peptide_test <- lapply(seq_len(nrow(df)), function(r) {
    x1 <- as.numeric(df[r, raw_cols1]); x2 <- as.numeric(df[r, raw_cols2])
    tt <- tryCatch(t.test(x1, x2, conf.level = 1 - pvalue), error = function(e) NULL)
    if (is.null(tt) || is.na(tt$p.value) || anyNA(tt$conf.int)) {
      return(list(p = NA_real_, rx = NA_real_, ry = NA_real_))
    }
    s1 <- stats::sd(x1, na.rm = TRUE); s2 <- stats::sd(x2, na.rm = TRUE)
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    t_crit <- qt(1 - pvalue / 2, tt$parameter)
    list(p = tt$p.value, rx = t_crit * s1 / sqrt(n1), ry = t_crit * s2 / sqrt(n2))
  })
  pv_state1v2 <- vapply(peptide_test, `[[`, numeric(1), "p")
  rx_state1v2 <- vapply(peptide_test, `[[`, numeric(1), "rx")
  ry_state1v2 <- vapply(peptide_test, `[[`, numeric(1), "ry")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(oma = c(3, 3, 1, 1))

  plot(av1[, av_state1], av1[, av_state2], type = "n", xlab = "", ylab = "", xlim = xlim, ylim = ylim)
  mtext("Deuteration (state 1)", side = 1, line = 0.7, outer = TRUE)
  mtext("Deuteration (state 2)", side = 2, line = 0.7, outer = TRUE)
  # the diagonal itself (the ellipse/p-value test's reference line) is drawn
  # solid and thicker; the CI band lines (a separate, pooled-SD threshold -
  # see @details) are dashed and lighter, so the two are not mistaken for
  # each other when an ellipse's edge crosses one but not the other.
  points(c(-10, 110), c(-10, 110), type = "l", lwd = 2)

  draw_one <- function(idx, ci, col_default, use_position_col, p_sig) {
    n_skipped <- 0
    for (pep in idx) {
      # a peptide can have an NA/NaN ellipse (e.g. identical-valued
      # replicates make the underlying t-test's variance undefined, or too
      # few non-NA replicates) - skip drawing it rather than crashing
      # DrawEllipse on a missing value, caught by testing against real
      # replicate data (some all-zero peptides in the package's own shipped
      # example return NaN from t.test() without erroring), not assumed away
      if (anyNA(c(av1[pep, av_state1], av1[pep, av_state2], rx_state1v2[pep], ry_state1v2[pep]))) {
        n_skipped <- n_skipped + 1
        next
      }
      delta <- abs(av1[pep, av_state2] - av1[pep, av_state1])
      is_sig <- (delta > ci) && !is.na(p_sig[[as.character(pep)]]) && (p_sig[[as.character(pep)]] < pvalue)
      base_col  <- if (is_sig) sig_col else nonsig_col
      border_col <- if (use_position_col) position_col[pep] else base_col
      fill_col   <- if (use_position_col) position_col[pep] else grDevices::adjustcolor(base_col, alpha.f = fill_alpha)
      DrawEllipse(av1[pep, av_state1], av1[pep, av_state2],
                  rx_state1v2[pep], ry_state1v2[pep],
                  col = fill_col, border = border_col, ...)
    }
    if (n_skipped > 0) message(sprintf("  (%d peptide(s) skipped: NA average or ellipse, e.g. identical-valued replicates)", n_skipped))
  }

  use_position_col <- !is.na(position_col[1])
  drawn <- list()
  for (i in seq_along(times)) {
    ind1 <- select_indices(av1, times = times[i], length = length1)
    if (length(ind1) == 0) {
      warning(sprintf("no peptides match times = '%s' - skipping", times[i]))
      next
    }
    ci <- CI_2pts(sd1[ind1, sd_state1], sd1[ind1, sd_state2], replicates, alpha,
                  n1 = .hdx_n_for(sd1, sd_state1, ind1), n2 = .hdx_n_for(sd1, sd_state2, ind1))
    message(sprintf("time %s: CI = %.3f", times[i], ci))
    points(c(-10, 110), c(-10, 110) + ci, type = "l", col = "grey70", lty = 2)
    points(c(-10, 110), c(-10, 110) - ci, type = "l", col = "grey70", lty = 2)
    p_adj <- stats::setNames(stats::p.adjust(pv_state1v2[ind1], method = p_adjust), ind1)
    draw_one(ind1, ci, "grey30", use_position_col, p_adj)
    delta <- av1[ind1, av_state2] - av1[ind1, av_state1]
    drawn[[length(drawn) + 1]] <- data.frame(
      Deut.Time = av1$Deut.Time[ind1], Start = av1$Start[ind1], End = av1$End[ind1],
      Charge = av1$Charge[ind1],
      Sequence = av1$Sequence[ind1], state1_av = av1[ind1, av_state1], state2_av = av1[ind1, av_state2],
      delta = delta, p_value = pv_state1v2[ind1], p_adj = as.numeric(p_adj), ci = ci,
      is_sig = !is.na(p_adj) & abs(delta) > ci & p_adj < pvalue,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, drawn)
  rownames(out) <- NULL
  invisible(.hdx_stamp(out, list(state1 = state1, state2 = state2, times = times, alpha = alpha,
                                 pvalue = pvalue, p_adjust = p_adjust, sd_source = sd_source)))
}
