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
#' @param palette either a single `RColorBrewer` qualitative palette name
#'   (default `"Set2"`) or an explicit vector of colors, one per state, for
#'   the peptide-count panel.
#' @param warn_col color for a shared-peptide bar below `low_threshold` in
#'   the overlap panel. Default `"firebrick"`.
#' @param ok_col color for a shared-peptide bar at or above `low_threshold`.
#'   Default `"grey55"`.
#' @param bg panel background color. Default `"grey95"`; use `"white"` for
#'   a plain white panel.
#' @param grid_col gridline color drawn over `bg`. Default `"white"`.
#' @param show_grid logical; draw the shaded background/gridlines at all.
#'   Default `TRUE`.
#' @param legend logical; draw the state legend on the count panel.
#'   Default `TRUE`.
#' @param cex_main text size for panel titles. Default `1`.
#' @param font_family font family passed to text-drawing calls. Default
#'   `""` (device default).
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
plot_coverage_diagnostics <- function(res, low_threshold = 0.5, palette = "Set2",
                                       warn_col = "firebrick", ok_col = "grey55",
                                       bg = "grey95", grid_col = "white", show_grid = TRUE,
                                       legend = TRUE, cex_main = 1, font_family = "") {
  cov <- res$coverage
  if (is.null(cov)) {
    stop("res$coverage not found - pass the list returned by compare_states_unmatched().")
  }
  regions <- unique(cov$region)
  states <- unique(cov$state)
  col <- .hdx_palette(length(states), palette)

  n_mat <- matrix(NA_real_, nrow = length(regions), ncol = length(states),
                   dimnames = list(regions, states))
  for (s in states) {
    sub <- cov[cov$state == s, ]
    n_mat[sub$region, s] <- sub$n_peptides
  }

  shared <- vapply(regions, function(r) cov$frac_shared_all_states[cov$region == r][1], numeric(1))
  bar_col <- ifelse(shared < low_threshold, warn_col, ok_col)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(family = font_family)
  par(mfrow = c(1, 2), mar = c(6, 4.5, 3.5, 1), oma = c(0, 0, 1, 0))

  ## generous headroom above the tallest bar when a legend is requested, so
  ## an in-panel top-right legend always lands in empty space regardless of
  ## how tall any individual bar happens to be
  ylim_n <- c(0, max(n_mat, na.rm = TRUE) * (if (legend) 1.35 else 1.08))
  graphics::barplot(t(n_mat), beside = TRUE, col = NA, border = NA, las = 2,
                     ylim = ylim_n, ylab = "Number of peptides", plot = TRUE)
  if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, ny = 5)
  graphics::barplot(t(n_mat), beside = TRUE, col = col, las = 2, ylim = ylim_n,
                     ylab = "Number of peptides", add = TRUE)
  .hdx_title("Peptides per region/state", cex_title = cex_main)
  if (legend) {
    graphics::legend("topright", legend = states, fill = col, bty = "n", cex = 0.85)
  }

  ylim_s <- c(0, 1)
  graphics::barplot(shared, names.arg = regions, col = NA, border = NA, las = 2,
                     ylim = ylim_s, ylab = "Fraction of peptide IDs shared by all states")
  if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, ny = 5)
  graphics::barplot(shared, names.arg = regions, col = bar_col, las = 2, ylim = ylim_s,
                     ylab = "Fraction of peptide IDs shared by all states", add = TRUE)
  graphics::abline(h = low_threshold, lty = 2, col = "grey30", lwd = 1.3)
  .hdx_title("Peptide-map overlap", cex_title = cex_main)

  invisible(NULL)
}

#' Plot the pooled distributions behind `compare_states_unmatched()`
#'
#' Draws what the comparison actually ran on: one panel per region, the
#' pooled `deut_value` values by state (one value per peptide, or per
#' replicate under `method = "classic"` - whatever `res$data` holds) as a
#' box, violin, or plain jittered-point strip (`plot_type`), with the
#' omnibus p-value and every pairwise result annotated on each panel. The
#' y-axis label follows `res$value` automatically (raw Da uptake by default -
#' see [compare_states_unmatched()]'s `value` argument). Works for any
#' number of states `compare_states_unmatched()` was run with, not only two.
#' This is the "show me the test" companion to
#' [plot_coverage_diagnostics()], which instead shows whether the test
#' should be trusted in the first place (peptide-map overlap) rather than
#' its result.
#'
#' @section On color and visibility:
#' Fills and points here are drawn as fully OPAQUE colors (box/violin fills
#' lightened by blending toward white in RGB space, not by alpha
#' transparency) deliberately - some graphics devices (notably Windows'
#' native `windows()` device, which RStudio's Plots pane can fall back to)
#' render semi-transparent fills or points as invisible or fully solid
#' instead of blending them, a real failure mode this avoids entirely by
#' never depending on device compositing for whether data is visible.
#'
#' @param res the list returned by [compare_states_unmatched()] (uses
#'   `res$data` and `res$summary`).
#' @param plot_type one of `"jitter"` (just the individual points, no
#'   summary shape - the safest default when you want to see every peptide,
#'   not only significant ones, at a glance), `"box"`, or `"violin"`
#'   (kernel density per state). Default `"jitter"`.
#' @param palette either a single `RColorBrewer` qualitative palette name
#'   (default `"Set2"`) or an explicit vector of colors, one per state.
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
#' @param show_points logical; overlay jittered individual data points on
#'   top of the box/violin shape. Always `TRUE` (and unused) when
#'   `plot_type = "jitter"`. Default `TRUE`.
#' @param point_cex size of the jittered points. Default `0.9`.
#' @param jitter_width horizontal jitter width for points, as a fraction of
#'   the box spacing. Default `0.15`.
#' @param tint how much lighter than the base state color the box/violin
#'   fill is, `0` (no lightening) to `1` (white) - see "On color and
#'   visibility". Default `0.55`.
#' @param notch logical; draw notched boxes (a visual approximation of a
#'   comparison interval around each median). Only applies to
#'   `plot_type = "box"`. Default `FALSE`.
#' @param bg panel background color. Default `"grey95"`; use `"white"` for
#'   a plain white panel.
#' @param grid_col gridline color drawn over `bg`. Default `"white"`.
#' @param show_grid logical; draw the shaded background/gridlines at all.
#'   Default `TRUE`.
#' @param font_family font family passed to text-drawing calls. Default
#'   `""` (device default).
#' @return draws a one-panel-per-region base-R plot; returns `NULL`
#'   invisibly.
#' @seealso [compare_states_unmatched()], [plot_coverage_diagnostics()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- compare_states_unmatched(file_nm, times = "1800.00s",
#'                                  regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' plot_compare_states_unmatched(res)
#' plot_compare_states_unmatched(res, plot_type = "violin",
#'                               state_labels = c(Unbound = "Apo", bound = "Complex"))
#' }
#' @export
plot_compare_states_unmatched <- function(res, plot_type = c("jitter", "box", "violin"),
                                            palette = "Set2", digits = 3, state_labels = NULL,
                                            cex_annot = 1, show_points = TRUE,
                                            point_cex = 0.9, jitter_width = 0.15, tint = 0.55,
                                            notch = FALSE, bg = "grey95", grid_col = "white",
                                            show_grid = TRUE, font_family = "") {
  plot_type <- match.arg(plot_type)
  data <- res$data
  summary_df <- res$summary
  if (is.null(data) || is.null(summary_df)) {
    stop("res$data/res$summary not found - pass the list returned by compare_states_unmatched().")
  }
  regions <- unique(data$region)
  states <- unique(data$Protein.State)
  col <- .hdx_palette(length(states), palette)
  names(col) <- states
  fill_col <- .hdx_tint(col, tint)
  names(fill_col) <- states

  value_label <- switch(if (is.null(res$value)) "raw" else res$value,
    raw = "Uptake (# Deut, Da)",
    percent = "% Deuterated",
    frac_deut = "Fraction deuterated (per residue)"
  )

  display <- function(s) {
    if (!is.null(state_labels) && s %in% names(state_labels)) unname(state_labels[s]) else s
  }
  display_labels <- vapply(states, display, character(1))

  ## enough top-margin lines for a bold region title plus one line per
  ## p-value (omnibus + one per pair), whichever region needs the most
  n_lines_max <- 3 + max(vapply(regions, function(r) sum(summary_df$region == r), integer(1)))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(family = font_family)
  par(mfrow = c(1, length(regions)), mar = c(7, 4.5, n_lines_max * 1.4 + 1, 1), oma = c(0, 0, 1, 0))

  for (r in regions) {
    sub <- data[data$region == r, ]
    sub$Protein.State <- factor(sub$Protein.State, levels = states, labels = display_labels)

    yr <- range(sub$deut_value, na.rm = TRUE)
    ylim <- yr + c(-1, 1) * diff(yr) * 0.08
    plot(NA, xlim = c(0.5, length(states) + 0.5), ylim = ylim, xaxt = "n",
         ylab = value_label, xlab = "")
    graphics::axis(1, at = seq_along(states), labels = display_labels, las = 2)
    if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, nx = 0, ny = 5)

    if (plot_type == "box") {
      graphics::boxplot(deut_value ~ Protein.State, data = sub, col = fill_col[states], border = "grey25",
                         staplewex = 0, outline = FALSE, add = TRUE, xaxt = "n", yaxt = "n")
      if (notch) {
        graphics::boxplot(deut_value ~ Protein.State, data = sub, notch = TRUE, col = NA, border = "grey15",
                           boxwex = 0.4, add = TRUE, xaxt = "n", yaxt = "n")
      }
    } else if (plot_type == "violin") {
      for (i in seq_along(states)) {
        yv <- sub$deut_value[sub$Protein.State == display_labels[i]]
        .hdx_violin(i, yv, col = fill_col[states[i]], border = "grey25")
      }
    }

    if (show_points || plot_type == "jitter") {
      for (i in seq_along(states)) {
        yv <- sub$deut_value[sub$Protein.State == display_labels[i]]
        xv <- i + stats::runif(length(yv), -jitter_width, jitter_width)
        graphics::points(xv, yv, pch = 16, cex = point_cex, col = col[states[i]])
      }
    }

    sr <- summary_df[summary_df$region == r, ]
    p_om <- sr$p_omnibus[1]
    body_lines <- sprintf("omnibus p = %s", signif(p_om, digits))
    for (i in seq_len(nrow(sr))) {
      body_lines <- c(body_lines, sprintf("%s vs %s: p_adj = %s", display(sr$state1[i]), display(sr$state2[i]),
                                            signif(sr$p_adj[i], digits)))
    }
    ## title sits farthest from the plot (highest line=), annotations stack
    ## below it in reading order, nearest line closest to the plot - so the
    ## region name always reads as the heading, not buried under its own
    ## p-value detail lines
    graphics::mtext(r, side = 3, line = length(body_lines) * 1.4 + 0.6, adj = 0, font = 2, cex = 1.1, col = "grey15")
    graphics::mtext(body_lines, side = 3, line = rev(seq(0.3, by = 1.4, length.out = length(body_lines))),
                     cex = cex_annot, adj = 0, col = "grey30")
  }

  invisible(NULL)
}
