#' Plot per-peptide deltas from `compare_states_tp()`, by state pair and region
#'
#' [compare_states_tp()] itself has no plotting companion - unlike
#' [compare_states_unmatched()] (paired with [plot_compare_states_unmatched()])
#' its peptide-by-pair result table was only ever inspected as a data.frame.
#' This draws one panel per region (or a single "all" panel if `regions`
#' wasn't used), with every matched peptide's `delta` shown as a box,
#' violin, or plain jittered-point strip - one position per state pair -
#' so a whole domain's behavior for a given comparison is visible at a
#' glance, for any number of states.
#'
#' @section Showing every peptide, not only the significant ones:
#' Every peptide tested is drawn, regardless of `is_sig` - significant ones
#' are simply colored differently (`sig_col` vs `nonsig_col`) rather than
#' filtered out. Looking only at a list of "significant peptides" hides
#' exactly the pattern [compare_states_unmatched()]'s
#' `frac_peptides_significant` diagnostic exists to catch: a region can
#' contain a real, consistent shift in most of its peptides while still
#' having only a few (or zero) individually clear the strict Hageman & Weis
#' gate, or vice versa - one dramatic outlier peptide surrounded by
#' otherwise flat ones. Seeing the whole distribution, not just the subset
#' that cleared a cutoff, is the point of this plot.
#'
#' @section On color and visibility:
#' Box/violin fills are opaque colors (lightened by blending toward white in
#' RGB space via `tint`, not by alpha transparency) deliberately - some
#' graphics devices render semi-transparent fills as invisible or fully
#' solid instead of blending them, so visibility here never depends on
#' device compositing.
#'
#' @param res the data.frame returned by [compare_states_tp()].
#' @param plot_type one of `"jitter"` (every peptide as a point, no summary
#'   shape), `"box"`, or `"violin"`. Default `"jitter"`.
#' @param sig_col color for peptides with `is_sig == TRUE`. Default
#'   `"firebrick"`.
#' @param nonsig_col color for peptides with `is_sig == FALSE`. Default
#'   `"grey55"`.
#' @param summary_col color for the box/violin summary shape (drawn once
#'   per pair, pooling significant and non-significant peptides together).
#'   Default `"steelblue"`.
#' @param tint how much lighter than `summary_col` the box/violin fill is,
#'   `0` (no lightening) to `1` (white). Default `0.6`.
#' @param point_cex size of the jittered points. Default `0.9`.
#' @param jitter_width horizontal jitter width for points, as a fraction of
#'   pair spacing. Default `0.15`.
#' @param show_zero logical; draw a dashed horizontal reference line at
#'   `delta = 0`. Default `TRUE`.
#' @param bg panel background color. Default `"grey95"`; use `"white"` for
#'   a plain white panel.
#' @param grid_col gridline color drawn over `bg`. Default `"white"`.
#' @param show_grid logical; draw the shaded background/gridlines at all.
#'   Default `TRUE`.
#' @param cex_annot text size for the "n significant / n tested" annotation
#'   under each region's title. Default `0.9`.
#' @param font_family font family passed to text-drawing calls. Default
#'   `""` (device default).
#' @return draws a one-panel-per-region base-R plot; returns `NULL`
#'   invisibly.
#' @seealso [compare_states_tp()], [plot_compare_states_unmatched()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' # the shipped example only has two states; compare_states_tp() needs 3+.
#' # See vignette("anova") for a worked example with a synthetic third state.
#' }
#' @export
plot_compare_states_tp <- function(res, plot_type = c("jitter", "box", "violin"),
                                    sig_col = "firebrick", nonsig_col = "grey55",
                                    summary_col = "steelblue", tint = 0.6,
                                    point_cex = 0.9, jitter_width = 0.15, show_zero = TRUE,
                                    bg = "grey95", grid_col = "white", show_grid = TRUE,
                                    cex_annot = 0.9, font_family = "") {
  plot_type <- match.arg(plot_type)
  if (!all(c("state1", "state2", "delta", "is_sig") %in% colnames(res))) {
    stop("res must be the data.frame returned by compare_states_tp().")
  }
  if (!"region" %in% colnames(res)) res$region <- "all"
  regions <- unique(res$region)
  res$pair <- paste(res$state1, "vs", res$state2)
  pairs <- unique(res$pair)
  fill_col <- .hdx_tint(summary_col, tint)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  graphics::par(family = font_family)
  par(mfrow = c(1, length(regions)), mar = c(8, 4.5, 3.5, 1), oma = c(0, 0, 1, 0))

  for (r in regions) {
    sub <- res[res$region == r, ]
    sub$pair <- factor(sub$pair, levels = pairs)

    yr <- range(sub$delta, na.rm = TRUE)
    if (show_zero) yr <- range(c(yr, 0))
    ylim <- yr + c(-1, 1) * diff(yr) * 0.1
    plot(NA, xlim = c(0.5, length(pairs) + 0.5), ylim = ylim, xaxt = "n",
         ylab = "delta", xlab = "")
    graphics::axis(1, at = seq_along(pairs), labels = pairs, las = 2)
    if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, nx = 0, ny = 5)
    if (show_zero) graphics::abline(h = 0, lty = 2, col = "grey30", lwd = 1.3)

    if (plot_type == "box") {
      graphics::boxplot(delta ~ pair, data = sub, col = fill_col, border = "grey25",
                         staplewex = 0, outline = FALSE, add = TRUE, xaxt = "n", yaxt = "n")
    } else if (plot_type == "violin") {
      for (i in seq_along(pairs)) {
        .hdx_violin(i, sub$delta[sub$pair == pairs[i]], col = fill_col, border = "grey25")
      }
    }

    for (i in seq_along(pairs)) {
      pd <- sub[sub$pair == pairs[i], ]
      pt_col <- ifelse(pd$is_sig, sig_col, nonsig_col)
      xv <- i + stats::runif(nrow(pd), -jitter_width, jitter_width)
      graphics::points(xv, pd$delta, pch = 16, cex = point_cex, col = pt_col)
    }

    n_sig <- sum(sub$is_sig, na.rm = TRUE)
    .hdx_title(r, sprintf("%d / %d peptide-pairs significant", n_sig, nrow(sub)), cex_sub = cex_annot)
    graphics::legend("topleft", legend = c("significant", "not significant"), col = c(sig_col, nonsig_col),
                      pch = 16, bty = "n", cex = 0.75)
  }

  invisible(NULL)
}
