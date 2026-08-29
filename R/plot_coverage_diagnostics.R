#' Plot the peptide-map coverage diagnostics from `compare_states_unmatched()`
#'
#' [compare_states_unmatched()] always returns a `coverage` diagnostic table
#' alongside its results, because its whole approach - pooling peptides
#' without requiring them to match between states - is only trustworthy
#' when the states being compared actually have similar peptide maps in
#' that region. This plots that table directly, two panels side by side:
#' peptide count per region/state (very different counts is itself a
#' warning sign), and the fraction of peptide identities shared by every
#' state being compared, per region (low values mean part of any observed
#' shift may reflect which peptides were identified, not biology - see
#' [compare_states_unmatched()]'s Details).
#'
#' @param res the list returned by [compare_states_unmatched()] (uses
#'   `res$coverage`; `res$summary` is not needed for this plot).
#' @param low_threshold fraction below which a region's shared-peptide bar
#'   is drawn in a warning color and a reference line is drawn. Default
#'   `0.5`.
#' @param col optional vector of colors, one per state, for the peptide
#'   count panel. Default `NULL` uses `RColorBrewer::brewer.pal(..., "Set2")`.
#' @return draws a two-panel base-R plot; returns `NULL` invisibly.
#' @seealso [compare_states_unmatched()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- compare_states_unmatched(file_nm, times = "1800.00s",
#'                                  regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' plot_coverage_diagnostics(res)
#' }
#' @export
plot_coverage_diagnostics <- function(res, low_threshold = 0.5, col = NULL) {
  cov <- res$coverage
  if (is.null(cov)) {
    stop("res$coverage not found - pass the list returned by compare_states_unmatched().")
  }
  regions <- unique(cov$region)
  states <- unique(cov$state)
  if (is.null(col)) {
    pal_n <- max(3, length(states))
    col <- RColorBrewer::brewer.pal(pal_n, "Set2")[seq_along(states)]
  }

  n_mat <- matrix(NA_real_, nrow = length(regions), ncol = length(states),
                   dimnames = list(regions, states))
  for (s in states) {
    sub <- cov[cov$state == s, ]
    n_mat[sub$region, s] <- sub$n_peptides
  }

  shared <- vapply(regions, function(r) cov$frac_shared_all_states[cov$region == r][1], numeric(1))
  bar_col <- ifelse(shared < low_threshold, "firebrick", "grey60")

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, 2), mar = c(6, 4, 3, 1), oma = c(0, 0, 1, 0))

  graphics::barplot(t(n_mat), beside = TRUE, col = col, las = 2,
                     ylab = "Number of peptides", main = "Peptides per region/state")
  graphics::legend("topright", legend = states, fill = col, bty = "n", cex = 0.8)

  graphics::barplot(shared, names.arg = regions, col = bar_col, ylim = c(0, 1), las = 2,
                     ylab = "Fraction of peptide IDs shared by all states",
                     main = "Peptide-map overlap")
  graphics::abline(h = low_threshold, lty = 2, col = "grey40")

  invisible(NULL)
}

#' Plot the pooled distributions behind `compare_states_unmatched()`
#'
#' Draws what the comparison actually ran on: one boxplot panel per region,
#' the pooled `frac_deut` values by state (one value per peptide, or per
#' replicate under `method = "classic"` - whatever `res$data` holds), with
#' the omnibus p-value and every pairwise result annotated on each panel.
#' This is the "show me the test" companion to
#' [plot_coverage_diagnostics()], which instead shows whether the test
#' should be trusted in the first place (peptide-map overlap) rather than
#' its result.
#'
#' @param res the list returned by [compare_states_unmatched()] (uses
#'   `res$data` and `res$summary`).
#' @param col optional vector of colors, one per state. Default `NULL` uses
#'   `RColorBrewer::brewer.pal(..., "Set2")`.
#' @param digits number of significant digits for the annotated p-values.
#'   Default `3`.
#' @param state_labels optional named character vector mapping raw state
#'   names (as they appear in the data, e.g. HDX Examiner's unedited
#'   `"Protein State 1"`) to display labels (e.g.
#'   `c("Protein State 1" = "Apo tubulin")`) - used on the x-axis and in
#'   the p-value annotations. States not named in this vector are shown
#'   as-is. Default `NULL` shows every state's raw name unchanged.
#' @param cex_annot text size for the p-value annotation lines. Default
#'   `1`.
#' @return draws a one-panel-per-region base-R boxplot; returns `NULL`
#'   invisibly.
#' @seealso [compare_states_unmatched()], [plot_coverage_diagnostics()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- compare_states_unmatched(file_nm, times = "1800.00s",
#'                                  regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' plot_compare_states_unmatched(res)
#' plot_compare_states_unmatched(res, state_labels = c(Unbound = "Apo", bound = "Complex"))
#' }
#' @export
plot_compare_states_unmatched <- function(res, col = NULL, digits = 3, state_labels = NULL,
                                            cex_annot = 1) {
  data <- res$data
  summary_df <- res$summary
  if (is.null(data) || is.null(summary_df)) {
    stop("res$data/res$summary not found - pass the list returned by compare_states_unmatched().")
  }
  regions <- unique(data$region)
  states <- unique(data$Protein.State)
  if (is.null(col)) {
    pal_n <- max(3, length(states))
    col <- RColorBrewer::brewer.pal(pal_n, "Set2")[seq_along(states)]
  }
  names(col) <- states

  display <- function(s) {
    if (!is.null(state_labels) && s %in% names(state_labels)) unname(state_labels[s]) else s
  }
  display_labels <- vapply(states, display, character(1))

  ## enough top-margin lines for a bold region title plus one line per
  ## p-value (omnibus + one per pair), whichever region needs the most
  n_lines_max <- 2 + max(vapply(regions, function(r) sum(summary_df$region == r), integer(1)))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, length(regions)), mar = c(7, 4.5, n_lines_max * 1.4 + 1, 1), oma = c(0, 0, 1, 0))

  for (r in regions) {
    sub <- data[data$region == r, ]
    sub$Protein.State <- factor(sub$Protein.State, levels = states, labels = display_labels)
    graphics::boxplot(frac_deut ~ Protein.State, data = sub, col = col[states],
                       ylab = "frac_deut", xlab = "", las = 2, cex.axis = 1)

    sr <- summary_df[summary_df$region == r, ]
    p_om <- sr$p_omnibus[1]
    body_lines <- sprintf("omnibus p = %s", signif(p_om, digits))
    for (i in seq_len(nrow(sr))) {
      body_lines <- c(body_lines, sprintf("%s vs %s: p_adj = %s", display(sr$state1[i]), display(sr$state2[i]),
                                            signif(sr$p_adj[i], digits)))
    }
    ## title (bold, larger) on its own line, farthest from the plot;
    ## p-value lines stacked below it, closest line nearest the plot box
    graphics::mtext(r, side = 3, line = length(body_lines) * 1.4 + 0.6, cex = 1.1, font = 2, adj = 0, col = "black")
    graphics::mtext(body_lines, side = 3, line = rev(seq(0.3, by = 1.4, length.out = length(body_lines))),
                     cex = cex_annot, adj = 0, col = "black")
  }

  invisible(NULL)
}
