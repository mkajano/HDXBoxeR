## Interactive HDX-onto-structure viewer. Optional-dependency function -
## bio3d and NGLVieweR are Suggests, not core dependencies - checked at call
## time with an install hint, the same pattern used in interactive_explorer.R
## for plotly/crosstalk/DT.

.hdx_require_structure <- function() {
  pkgs <- c("bio3d", "NGLVieweR")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Package(s) %s required for plot_hdx_structure() but not installed. Install with install.packages(c(%s)).",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf("'%s'", missing_pkgs), collapse = ", ")
    ), call. = FALSE)
  }
}

#' Map residue-level HDX data onto a 3D structure
#'
#' Colors a PDB structure by per-residue HDX difference (`state2 - state1`
#' at one timepoint), the structural counterpart to [heat_map_tp()]-style
#' 2D heat maps: same red/blue diverging palette
#' ([color_ranges_Blue_Red_heat_map()]), but rendered on the actual fold
#' instead of a sequence axis. This is an R-native alternative to the
#' package's existing [pymol_script_average_residue()] family, which writes
#' a text script that must be run inside a separately-installed copy of
#' PyMOL against a structure the user opens themselves; this function
#' instead renders directly as an interactive `NGLVieweR`/`NGL.js`
#' `htmlwidget` - viewable in the RStudio viewer, embeddable in an R
#' Markdown document, or saved standalone with `htmlwidgets::saveWidget()`.
#'
#' @param res_df residue-level data.frame - the output of
#'   [residue_level_average()] or [residue_deconvolution()] (one row per
#'   residue and `Deut.Time`, one column per Protein.State/replicate).
#' @param pdb a local PDB file path, or a 4-character PDB ID (fetched by
#'   `NGLVieweR`/`bio3d` directly, e.g. `"1HEL"`) - whatever
#'   `NGLVieweR::NGLVieweR()` itself accepts as `data`.
#' @param state1 First protein state, matched the same way as in
#'   [plot_fissure_vent()] (substring match against `res_df`'s column
#'   names).
#' @param state2 Second protein state, same matching rule.
#' @param times a single `Deut.Time` value to color by.
#' @param chain optional PDB chain identifier (e.g. `"A"`) to restrict
#'   coloring to. `NULL` (default) colors matching residue numbers in every
#'   chain.
#' @param resno_offset integer added to `res_df`'s residue numbers before
#'   matching them against the PDB's residue numbering. Default `0`. HDX
#'   peptide numbering and PDB numbering frequently disagree (constructs
#'   with tags, truncated crystallized regions, etc.) - this is a constant
#'   shift, not a general remapping; a construct that doesn't differ from
#'   the PDB by a single constant offset needs its `res_df$residue` column
#'   remapped by hand before calling this function.
#' @param ranges bin edges for the diverging color scale, passed to
#'   [color_ranges_Blue_Red_heat_map()] - must have odd length and be
#'   symmetric around (and include) `0`, e.g.
#'   `c(-Inf, -20, -10, 0, 10, 20, Inf)`. Default `NULL` auto-scales a
#'   7-bin symmetric range from the largest `|delta|` present, in whatever
#'   units `res_df` is already in (Da, percent, or fraction - unchanged
#'   from `res_df`, not re-derived here).
#' @param col color for residues with no HDX data at this
#'   timepoint (`NA` delta, or simply absent from `res_df`). Default
#'   `"grey80"`.
#' @param colors_initial deprecated, use `col`.
#' @param background background color for the 3D viewer. Default `"white"`.
#' @return an `NGLVieweR` `htmlwidget`.
#' @seealso [residue_level_average()], [residue_deconvolution()],
#'   [color_ranges_Blue_Red_heat_map()], [pymol_script_average_residue()]
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE) &&
#'     requireNamespace("NGLVieweR", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm, percent = TRUE)
#' res <- residue_level_average(a, protein_length = 225)
#' pdb_demo <- system.file("examples", "1hel.pdb", package = "bio3d")
#' plot_hdx_structure(res, pdb = pdb_demo, state1 = "Unbound", state2 = "bound",
#'                     times = "1800.00s")
#' }
#' @export
plot_hdx_structure <- function(res_df, pdb, state1, state2, times,
                                chain = NULL, resno_offset = 0,
                                ranges = NULL, col = "grey80",
                                background = "white", colors_initial) {
  if (!missing(colors_initial)) { .hdx_deprecated("colors_initial", "col"); col <- colors_initial }
  .hdx_require_structure()

  sub_df <- res_df[res_df$Deut.Time == times, ]
  if (nrow(sub_df) == 0) {
    stop(sprintf("no rows with Deut.Time == '%s'", times))
  }

  # residue_level_average()/residue_deconvolution() return one column per
  # REPLICATE (same naming as the raw X..Deut/Deut.. columns), not one
  # pre-averaged column per state - average whatever matches here rather
  # than requiring exactly one column, or this breaks on any normal
  # multi-replicate res_df.
  resolve_state_mean <- function(state, cols, mat) {
    ## exact state match (an unambiguous abbreviation is accepted), so "WT"
    ## does not also pick up "WT_mut"
    hits <- .hdx_state_cols(state, cols)
    if (length(hits) == 1) return(mat[[hits]])
    rowMeans(as.matrix(mat[, hits, drop = FALSE]), na.rm = TRUE)
  }
  state1_val <- resolve_state_mean(state1, colnames(sub_df), sub_df)
  state2_val <- resolve_state_mean(state2, colnames(sub_df), sub_df)

  delta <- state2_val - state1_val
  pdb_resno <- as.integer(sub_df$residue) + resno_offset

  if (is.null(ranges)) {
    m <- max(abs(delta), na.rm = TRUE)
    if (!is.finite(m) || m == 0) {
      stop("cannot auto-scale `ranges` (delta is all NA/zero) - supply `ranges` explicitly.")
    }
    step <- m / 3
    ranges <- c(-Inf, -3 * step, -2 * step, -step, 0, step, 2 * step, 3 * step, Inf)
  }

  # validate pdb/chain/numbering with bio3d before rendering, so a mismatch
  # (wrong chain, or HDX numbering that doesn't line up with the PDB) fails
  # loudly instead of silently producing an all-grey structure - bio3d only
  # needs to be able to read a LOCAL file; a 4-character PDB ID handed
  # straight to NGLVieweR (which fetches it itself) skips this check.
  pdb_obj <- if (file.exists(pdb)) tryCatch(bio3d::read.pdb(pdb, verbose = FALSE), error = function(e) NULL) else NULL
  if (!is.null(pdb_obj)) {
    ca <- pdb_obj$atom[pdb_obj$calpha, ]
    if (!is.null(chain) && !(chain %in% ca$chain)) {
      stop(sprintf("chain '%s' not found in `pdb` - chains present: %s",
                    chain, paste(unique(ca$chain), collapse = ", ")))
    }
    struct_resno <- if (is.null(chain)) ca$resno else ca$resno[ca$chain == chain]
    if (sum(pdb_resno %in% struct_resno, na.rm = TRUE) == 0) {
      warning(sprintf(
        "none of the %d residue numbers (after resno_offset = %d) were found in `pdb`%s - check numbering with `resno_offset`.",
        length(pdb_resno), resno_offset,
        if (!is.null(chain)) paste0(" chain '", chain, "'") else ""), call. = FALSE)
    }
  }

  colors_all <- color_ranges_Blue_Red_heat_map(ranges, col)
  bin_idx <- findInterval(delta, ranges)
  bin_idx[is.na(bin_idx) | bin_idx < 1 | bin_idx > length(ranges) - 1] <- 0

  ngl <- NGLVieweR::NGLVieweR(pdb) %>%
    NGLVieweR::stageParameters(backgroundColor = background) %>%
    NGLVieweR::addRepresentation("cartoon", param = list(name = "base", color = colors_all[1]))

  for (b in seq_len(length(ranges) - 1)) {
    resi_b <- pdb_resno[which(bin_idx == b)]
    resi_b <- resi_b[!is.na(resi_b)]
    if (length(resi_b) == 0) next
    sele <- if (is.null(chain)) paste(resi_b, collapse = ", ") else paste0(resi_b, ":", chain, collapse = ", ")
    ngl <- ngl %>% NGLVieweR::addRepresentation("cartoon", param = list(sele = sele, color = colors_all[b + 1]))
  }
  ngl
}

#' Legend for [plot_hdx_structure()]'s color scale
#'
#' Thin wrapper around [pallette_ll()] - draws the same red/blue bins used
#' by [plot_hdx_structure()] as a base-R legend, since the color scale isn't
#' otherwise labeled inside the 3D `NGLVieweR` widget itself.
#'
#' @param ranges the same `ranges` vector passed to [plot_hdx_structure()]
#'   (or, if that call used the default `NULL`, the auto-scaled vector -
#'   see [plot_hdx_structure()]'s `@details` on reconstructing it).
#' @param col the same `col` passed to
#'   [plot_hdx_structure()]. Default `"grey80"`.
#' @param colors_initial deprecated, use `col`.
#' @param unit_label axis label for the bins, e.g. `"Da"` or `"%"`. Default
#'   `""`.
#' @return draws a base-R legend plot; returns `NULL` invisibly.
#' @seealso [plot_hdx_structure()]
#' @examples
#' plot_hdx_structure_legend(ranges = c(-Inf, -20, -10, 0, 10, 20, Inf), unit_label = "Da")
#' @export
plot_hdx_structure_legend <- function(ranges, col = "grey80", unit_label = "", colors_initial) {
  if (!missing(colors_initial)) { .hdx_deprecated("colors_initial", "col"); col <- colors_initial }
  colors_all <- color_ranges_Blue_Red_heat_map(ranges, col)
  lab <- c("no data")
  for (i in seq_len(length(ranges) - 1)) {
    lab <- c(lab, sprintf("%s:%s %s", round(ranges[i], 1), round(ranges[i + 1], 1), unit_label))
  }
  pallette_ll(colors_all, lab)
  invisible(NULL)
}
