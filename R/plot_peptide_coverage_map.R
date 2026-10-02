## internal helper for plot_peptide_coverage_map() - not exported

## Greedy first-fit row packing: assigns each peptide (given as parallel
## Start/End vectors, any order) to the lowest-numbered row whose
## previously-placed peptide ends at least `gap` residues before this one
## starts, or a new row if none qualifies. Standard interval-graph coloring
## heuristic - not guaranteed minimum rows, but always visually
## non-overlapping and fast enough for typical peptide counts (hundreds).
.pack_peptide_rows <- function(starts, ends, gap = 1) {
  ord <- order(starts, ends)
  row_end <- numeric(0)
  row_of <- integer(length(starts))
  for (i in ord) {
    placed <- FALSE
    for (r in seq_along(row_end)) {
      if (row_end[r] + gap < starts[i]) {
        row_of[i] <- r
        row_end[r] <- ends[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      row_end <- c(row_end, ends[i])
      row_of[i] <- length(row_end)
    }
  }
  row_of
}

#' Plot peptide coverage map with per-residue redundancy
#'
#' The classic HDX-MS QC figure that this package didn't have a plotting
#' function for yet: [general_info()] already reports peptide coverage %,
#' peptide count, and mean redundancy as single numbers per state, but a
#' number can't show *where* a coverage gap is or *which* residues are
#' over/under-represented - only a map can. For each state, draws two
#' aligned panels sharing one residue axis: a per-residue redundancy strip
#' (how many peptides cover each residue - built from [coverage_residue()])
#' above a peptide map (one block per unique peptide, positioned at its
#' Start-End span, stacked into just enough rows to avoid overlap).
#'
#' Peptide identity here is `Start`/`End`/`Sequence`, deduplicated across
#' every `Deut.Time` and replicate for that state (a peptide's identity
#' doesn't change across timepoints, only whether it was successfully
#' quantified there) - so this shows what was *identified* for a state, the
#' same question [compare_states_unmatched()]'s `coverage`/
#' [plot_coverage_diagnostics()] ask numerically. Comparing two states'
#' panels side by side (same x-axis) is a direct visual check of exactly
#' the confound that function's `frac_shared_all_states` warns about: a gap
#' in one state's map that isn't a gap in another's marks where an unmatched
#' comparison could be measuring peptide-map differences, not biology.
#'
#' @param filepath filepath to the input file (read via
#'   [all_data_long_input()]).
#' @param states character vector of states to show, one redundancy+map
#'   panel pair per state, stacked top to bottom on a shared residue axis.
#'   Default `NULL` uses every state found in `filepath`.
#' @param palette either a single `RColorBrewer` qualitative palette name
#'   (default `"Set2"`) or an explicit vector of colors, one per state (or
#'   more; the first `length(states)` are used).
#' @param bg panel background color (the shaded area behind the gridlines
#'   and data). Default `"grey95"`; use `"white"` for a plain white panel.
#' @param grid_col gridline color drawn over `bg`. Default `"white"`.
#' @param show_grid logical; draw the shaded background/gridlines at all.
#'   Default `TRUE`; set `FALSE` for a plain white panel with no grid.
#' @param row_gap minimum residue gap enforced between two peptides placed
#'   in the same map row - visual spacing only, does not affect the
#'   redundancy strip. Default `1`.
#' @param block_height height of each peptide block as a fraction of one
#'   row (`0`-`1`). Default `0.5` - leaves visible spacing between rows.
#' @param block_border border color for peptide blocks, or `NA` for no
#'   border. Default `NA`.
#' @param cex_main text size for each state's title line (state name,
#'   peptide count, percent sequence coverage). Default `1`.
#' @param font_family font family passed to text-drawing calls (titles,
#'   axis labels). Default `""` (device default).
#' @return draws a two-panel-per-state base-R plot; returns `NULL`
#'   invisibly.
#' @seealso [general_info()], [coverage_residue()],
#'   [compare_states_unmatched()], [plot_coverage_diagnostics()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' plot_peptide_coverage_map(file_nm)
#' plot_peptide_coverage_map(file_nm, states = "Unbound", palette = "steelblue", bg = "white")
#' }
#' @export
plot_peptide_coverage_map <- function(filepath, states = NULL, palette = "Set2",
                                       bg = "grey95", grid_col = "white", show_grid = TRUE,
                                       row_gap = 1, block_height = 0.5, block_border = NA,
                                       cex_main = 1, font_family = "") {
  b <- all_data_long_input(filepath)
  if (is.null(states)) states <- unique(b$Protein.State)
  if (!all(states %in% unique(b$Protein.State))) {
    stop(sprintf("state(s) not found: %s", paste(setdiff(states, unique(b$Protein.State)), collapse = ", ")))
  }
  col <- .hdx_palette(length(states), palette)
  names(col) <- states

  peps <- lapply(states, function(s) {
    d <- unique(b[b$Protein.State == s, c("Start", "End", "Sequence")])
    d[order(d$Start, d$End), ]
  })
  names(peps) <- states

  xmax <- max(vapply(peps, function(d) max(d$End), numeric(1)))

  redund <- lapply(peps, function(d) {
    cov <- coverage_residue(d, start_col = 1, end_col = 2)
    if (length(cov) < xmax) cov <- c(cov, rep(0, xmax - length(cov)))
    cov
  })

  oldpar <- par(no.readonly = TRUE)
  on.exit({ par(oldpar); graphics::layout(1) })
  graphics::par(family = font_family)
  graphics::layout(matrix(seq_len(2 * length(states)), ncol = 1), heights = rep(c(1, 2.2), length(states)))

  for (i in seq_along(states)) {
    s <- states[i]
    cov <- redund[[s]]
    d <- peps[[s]]
    row_of <- .pack_peptide_rows(d$Start, d$End, gap = row_gap)
    pct_cov <- round(100 * sum(cov > 0) / length(cov), 1)
    last <- i == length(states)

    ## redundancy strip: filled area, not a needle plot - reads more like a
    ## coverage-depth track (genome-browser convention) at a glance
    par(mar = c(0.3, 4.5, 2.6, 1))
    ylim_d <- c(0, max(cov, 1) * 1.05)
    plot(NA, xlim = c(1, xmax), ylim = ylim_d, xaxs = "i", yaxs = "i", xaxt = "n", ylab = "depth", xlab = "")
    if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, ny = 3)
    graphics::polygon(c(1, seq_along(cov), xmax), c(0, cov, 0), col = col[s], border = NA)
    .hdx_title(sprintf("%s", s), sprintf("%d peptides · %s%% sequence coverage", nrow(d), pct_cov),
               cex_title = cex_main)

    ## peptide map: solid blocks (rect), not thin line segments - reads as
    ## a stacked coverage track rather than a scatter of dashes
    par(mar = c(if (last) 3.5 else 0.3, 4.5, 0.3, 1))
    ylim_m <- c(0.5, max(row_of) + 0.5)
    plot(NA, xlim = c(1, xmax), ylim = ylim_m, xaxs = "i", yaxt = "n",
         ylab = "peptides", xlab = if (last) "Residue" else "",
         xaxt = if (last) "s" else "n")
    if (show_grid) .hdx_panel(par("usr")[1:2], par("usr")[3:4], bg = bg, grid_col = grid_col, ny = 0)
    h <- block_height / 2
    graphics::rect(d$Start, row_of - h, d$End, row_of + h,
                    col = col[s], border = block_border)
  }

  invisible(NULL)
}
