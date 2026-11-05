#' Compare two states over a whole time course
#'
#' Combines the labeling times of each peptide into whole-time-course
#' statistics, so a peptide is judged on all its time points together
#' rather than one time point at a time.
#'
#' For each peptide (present at every requested time in both states):
#' \describe{
#'   \item{Summed difference}{`summed_delta` = sum over time points of
#'     (mean state2 - mean state1), in Da (or percent). Its standard error
#'     combines the per-time standard errors, `p_value` is a Welch-type
#'     t-test of the sum (Welch-Satterthwaite degrees of freedom), and
#'     `ci_sum` is a global critical interval for summed differences: the
#'     per-time pooled Hageman & Weis standard errors combined over time
#'     points, times `|qt(alpha/2, df)|`. `sig_summed`
#'     uses the package-wide rule: (adjusted) p < `pvalue` AND
#'     |summed_delta| > `ci_sum`.}
#'   \item{Per time point}{the usual peptide rule at each time (Welch p,
#'     adjusted across peptides within that time with `p_adjust`, below
#'     `pvalue` AND |delta| > that time's pooled CI); `n_timepoints_sig` is
#'     how many time points pass, and `sig_min_timepoints` whether that is
#'     at least `min_timepoints`.}
#' }
#'
#' @param filepath input file or table, as for [output_tp()].
#' @param state1,state2 the two protein states (differences are
#'   `state2 - state1`).
#' @param times time points to combine. Default: all common time points.
#' @param alpha alpha of the critical intervals. Default `0.01`.
#' @param pvalue p-value cutoff. Default `0.01`.
#' @param p_adjust multiple-testing correction passed to [stats::p.adjust()]:
#'   `"none"` (default) or e.g. `"BH"`; across peptides for the summed test,
#'   and across peptides within each time for the per-time calls.
#' @param min_timepoints number of time points at which a peptide must pass
#'   the per-time rule for `sig_min_timepoints`. Default `2`.
#' @param percent value scale as in [output_tp()]: `FALSE` (Da, default),
#'   `TRUE` or `"FD"`.
#' @param replicates optional maximum number of replicates per state.
#' @return data frame with one row per peptide: `Start`, `End`, `Sequence`,
#'   `Charge`, `n_times`, `summed_delta`, `se`, `df`, `p_value`, `p_adj`,
#'   `ci_sum`, `sig_summed`, `n_timepoints_sig`, `sig_min_timepoints`, and
#'   one `sig_<time>` column per time point.
#' @seealso [robot_indexes()], [woods_CI_plot()], [compare_regions()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' tc <- compare_timecourse(file_nm, "Unbound", "bound")
#' head(tc[order(tc$p_value), ])
#' @export
compare_timecourse <- function(filepath, state1, state2, times, alpha = 0.01, pvalue = 0.01,
                               p_adjust = "none", min_timepoints = 2, percent = FALSE,
                               replicates = NULL) {
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  raw <- .hdx_raw(filepath)
  a <- suppressMessages(if (missing(times)) {
    output_tp(raw, states = c(state1, state2), replicates = replicates, percent = percent)
  } else {
    output_tp(raw, states = c(state1, state2), times = times, replicates = replicates, percent = percent)
  })
  tlev <- unique(a$Deut.Time)
  m1 <- .hdx_block_matrix(a, .hdx_state_cols(state1, colnames(a)))
  m2 <- .hdx_block_matrix(a, .hdx_state_cols(state2, colnames(a)))
  a$mu1 <- .hdx_row_mean(m1); a$mu2 <- .hdx_row_mean(m2)
  a$s1 <- .hdx_row_sd(m1); a$s2 <- .hdx_row_sd(m2)
  a$n1 <- .hdx_row_n(m1); a$n2 <- .hdx_row_n(m2)
  a$delta <- a$mu2 - a$mu1
  a$p <- .hdx_row_welch(m1, m2)
  a$v <- a$s1^2 / a$n1 + a$s2^2 / a$n2
  a$pep <- paste(a$Start, a$End, a$Sequence, a$Charge, sep = "\r")

  ## per-time pooled SE and CI (Hageman & Weis) and per-time calls
  se_t <- ci_t <- df_t <- stats::setNames(numeric(length(tlev)), tlev)
  a$sig <- FALSE
  for (t in tlev) {
    sel <- a$Deut.Time == t
    ok <- sel & !is.na(a$s1) & !is.na(a$s2)
    sp1 <- .hdx_pooled_sd(a$s1[ok], a$n1[ok]); sp2 <- .hdx_pooled_sd(a$s2[ok], a$n2[ok])
    N1 <- mean(a$n1[ok]); N2 <- mean(a$n2[ok])
    se_t[[t]] <- sqrt(sp1^2 / N1 + sp2^2 / N2)
    df_t[[t]] <- N1 + N2 - 2
    ci_t[[t]] <- se_t[[t]] * abs(stats::qt(alpha/2, df_t[[t]]))
    padj <- stats::p.adjust(a$p[sel], method = p_adjust)
    a$sig[sel] <- !is.na(padj) & padj < pvalue & !is.na(a$delta[sel]) & abs(a$delta[sel]) > ci_t[[t]]
  }

  peps <- unique(a$pep)
  rows <- lapply(peps, function(pp) {
    d <- a[a$pep == pp, , drop = FALSE]
    complete <- all(tlev %in% d$Deut.Time) && !anyNA(d$v)
    sd_sum <- sum(d$delta)
    se <- sqrt(sum(d$v))
    dfw <- sum(d$v)^2 / sum((d$s1^2 / d$n1)^2 / (d$n1 - 1) + (d$s2^2 / d$n2)^2 / (d$n2 - 1))
    pval <- if (complete && se > 0) 2 * stats::pt(-abs(sd_sum / se), dfw) else NA_real_
    sig_cols <- stats::setNames(as.list(tlev %in% d$Deut.Time[d$sig]), paste0("sig_", tlev))
    data.frame(Start = d$Start[1], End = d$End[1], Sequence = d$Sequence[1], Charge = d$Charge[1],
               n_times = nrow(d),
               summed_delta = if (complete) sd_sum else NA_real_,
               se = if (complete) se else NA_real_, df = if (complete) dfw else NA_real_,
               p_value = pval, n_timepoints_sig = sum(d$sig), sig_cols,
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  ci_sum <- sqrt(sum(se_t^2)) * abs(stats::qt(alpha/2, sum(df_t)))
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust)
  out$ci_sum <- ci_sum
  out$sig_summed <- !is.na(out$p_adj) & out$p_adj < pvalue & abs(out$summed_delta) > ci_sum
  out$sig_min_timepoints <- out$n_timepoints_sig >= min_timepoints
  first <- c("Start", "End", "Sequence", "Charge", "n_times", "summed_delta", "se", "df",
             "p_value", "p_adj", "ci_sum", "sig_summed", "n_timepoints_sig", "sig_min_timepoints")
  out <- out[order(out$Start, out$End, out$Charge), c(first, setdiff(colnames(out), first))]
  rownames(out) <- NULL
  attr(out, "times") <- tlev
  .hdx_stamp(out, list(state1 = state1, state2 = state2, times = tlev, alpha = alpha, pvalue = pvalue,
                       p_adjust = p_adjust, min_timepoints = min_timepoints, percent = format(percent)))
}
