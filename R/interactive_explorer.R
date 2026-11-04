## Interactive volcano / fissure_vent / Woods plots, linked to a shared
## significant-peptide table via crosstalk. These are optional-dependency
## functions (plotly, crosstalk, DT) - the base package install does not
## require them; each exported function below checks for them at call time
## and errors with an install hint if missing, the standard CRAN pattern
## for Suggests-only functionality.

.hdx_require_interactive <- function(need_dt = FALSE, need_htmltools = FALSE) {
  pkgs <- c("plotly", "crosstalk")
  if (need_dt) pkgs <- c(pkgs, "DT")
  if (need_htmltools) pkgs <- c(pkgs, "htmltools")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf(
      "Package(s) %s required for interactive plots but not installed. Install with install.packages(c(%s)).",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf("'%s'", missing_pkgs), collapse = ", ")
    ), call. = FALSE)
  }
}

## Builds the one per-peptide table every interactive function in this file
## is a view onto: average deuteration in both states, the per-peptide
## Welch delta/p-value, and matching ellipse radii - the exact same
## statistics as plot_fissure_vent() (state1 vs state2, that peptide's own
## raw replicates, conf.level = 1 - pvalue so "ellipse touches zero" is
## exactly "p >= pvalue" - see vignette("fissure_vent")), computed
## independently here rather than by calling plot_fissure_vent() (which
## draws a base-R plot as a side effect and isn't a data-returning function).
.hdx_comparison_table <- function(df, state1, state2, times, replicates, alpha, pvalue) {
  sub_df <- df[df$Deut.Time == times, ]
  if (nrow(sub_df) == 0) {
    stop(sprintf("no rows with Deut.Time == '%s'", times))
  }

  av1 <- ave_timepoint(sub_df, replicates)
  nm1 <- sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", sub("^av_", "", colnames(av1)[7:ncol(av1)]))
  if (missing(state1) || is.null(state1)) state1 <- nm1[1]
  if (missing(state2) || is.null(state2)) state2 <- nm1[2]

  resolve_col <- function(state, cols, known_prefixes = character(0)) {
    state <- .state_column_prefix(state)
    root <- cols
    for (p in known_prefixes) root <- sub(paste0("^", p), "", root)
    hits <- cols[startsWith(root, state)]
    if (length(hits) != 1) {
      stop(sprintf("State '%s' must match exactly one column among: %s",
                    state, paste(cols, collapse = ", ")))
    }
    hits
  }
  av_state1 <- resolve_col(state1, colnames(av1), "av_")
  av_state2 <- resolve_col(state2, colnames(av1), "av_")

  ## raw per-state replicate SD - same pooled CI as plot_fissure_vent()'s
  ## default and the rest of the package (CI_tp(), diff_peptides_prep())
  sd1 <- sd_timepoint(sub_df, replicates)
  sd_prefixes <- c("sd_ci_", "sd_")
  sd_state1 <- resolve_col(state1, colnames(sd1), sd_prefixes)
  sd_state2 <- resolve_col(state2, colnames(sd1), sd_prefixes)
  ci <- CI_2pts(stats::na.omit(sd1[[sd_state1]]), stats::na.omit(sd1[[sd_state2]]), replicates, alpha)

  resolve_raw_cols <- function(state, cols, replicates) {
    state <- .state_column_prefix(state)
    hits <- cols[startsWith(cols, state)]
    if (length(hits) != replicates) {
      stop(sprintf("Expected %d raw replicate columns matching state '%s', found %d: %s",
                    replicates, state, length(hits), paste(hits, collapse = ", ")))
    }
    hits
  }
  raw_cols1 <- resolve_raw_cols(state1, colnames(sub_df), replicates)
  raw_cols2 <- resolve_raw_cols(state2, colnames(sub_df), replicates)

  peptide_test <- lapply(seq_len(nrow(sub_df)), function(r) {
    x1 <- as.numeric(sub_df[r, raw_cols1]); x2 <- as.numeric(sub_df[r, raw_cols2])
    tt <- tryCatch(t.test(x1, x2, conf.level = 1 - pvalue), error = function(e) NULL)
    if (is.null(tt) || is.na(tt$p.value) || anyNA(tt$conf.int)) {
      return(list(p = NA_real_, rx = NA_real_, ry = NA_real_))
    }
    s1 <- stats::sd(x1, na.rm = TRUE); s2 <- stats::sd(x2, na.rm = TRUE)
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    t_crit <- qt(1 - pvalue / 2, tt$parameter)
    list(p = tt$p.value, rx = t_crit * s1 / sqrt(n1), ry = t_crit * s2 / sqrt(n2))
  })

  tbl <- data.frame(
    peptide_id = paste(av1$Start, av1$End, av1$Sequence, sep = "_"),
    Start = av1$Start, End = av1$End, Sequence = av1$Sequence,
    state1_val = av1[[av_state1]], state2_val = av1[[av_state2]],
    delta = av1[[av_state2]] - av1[[av_state1]],
    p_value = vapply(peptide_test, `[[`, numeric(1), "p"),
    rx = vapply(peptide_test, `[[`, numeric(1), "rx"),
    ry = vapply(peptide_test, `[[`, numeric(1), "ry"),
    ci = ci,
    stringsAsFactors = FALSE
  )
  tbl$is_sig <- !is.na(tbl$p_value) & (abs(tbl$delta) > tbl$ci) & (tbl$p_value < pvalue)
  tbl$hover <- sprintf("%s\nStart-End: %d-%d\n%s: %.2f  %s: %.2f\ndelta: %.2f   p = %.3g%s",
                        tbl$Sequence, tbl$Start, tbl$End, state1, tbl$state1_val,
                        state2, tbl$state2_val, tbl$delta, tbl$p_value,
                        ifelse(tbl$is_sig, "  [significant]", ""))
  attr(tbl, "state1") <- state1
  attr(tbl, "state2") <- state2
  attr(tbl, "ci") <- ci
  attr(tbl, "pvalue") <- pvalue
  tbl
}

.hdx_shared <- function(df, state1, state2, times, replicates, alpha, pvalue, shared_data, group_prefix) {
  if (!is.null(shared_data)) return(shared_data)
  tbl <- .hdx_comparison_table(df, state1, state2, times, replicates, alpha, pvalue)
  crosstalk::SharedData$new(tbl, key = ~peptide_id,
                             group = paste0(group_prefix, "_", paste(sample(letters, 8), collapse = "")))
}

#' Interactive volcano plot, linked to the fissure_vent/Woods/table views
#'
#' Plotly version of [plots_vol_tp()]'s delta-vs-p-value display: x = state2
#' - state1 (Da or the units of `df`), y = -log10(p-value), with the same
#' hybrid criterion (Hageman & Weis pooled CI + per-peptide Welch p-value)
#' used everywhere else in this package. Hover shows the peptide sequence,
#' range, both states' values, delta, and p-value.
#'
#' @param df output of [output_tp()].
#' @param state1 First protein state, matched as in [plot_fissure_vent()].
#'   Defaults to the first state found in `df`.
#' @param state2 Second protein state. Defaults to the second state found.
#' @param times a single Deut.Time value to plot.
#' @param replicates number of replicates. Default `3`.
#' @param alpha alpha for the pooled critical interval. Default `0.01`. May
#'   be a vector (paired element-wise with `pvalue`, recycled to a common
#'   length) - with more than one combination, the plot gets a dropdown
#'   (native `plotly` `updatemenus`, no R session needed) that swaps between
#'   the pre-computed colorings/thresholds for each `alpha`/`pvalue` pair.
#'   Peptide *positions* never move (delta and p-value don't depend on
#'   either parameter) - only which points/lines count as significant.
#' @param pvalue per-peptide significance cutoff. Default `0.01`. See
#'   `alpha` for vector/dropdown behavior.
#' @param shared_data optional `crosstalk::SharedData` object (as built
#'   internally by [interactive_fissure_explorer()]) so this plot's
#'   selection is linked to other plots sharing the same object. If `NULL`
#'   (default), a new one is built and this plot works standalone. The
#'   linked *key* (which peptide is which) always follows `shared_data`;
#'   only the significance coloring/dropdown below is local to this plot.
#' @return a `plotly` htmlwidget.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' interactive_volcano(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' interactive_volcano(a, state1 = "Unbound", state2 = "bound", times = "1800.00s",
#'                      alpha = c(0.01, 0.05), pvalue = c(0.01, 0.05))
#' }
#' @export
interactive_volcano <- function(df, state1, state2, times, replicates = 3,
                                 alpha = 0.01, pvalue = 0.01, shared_data = NULL) {
  .hdx_require_interactive()
  n_combo <- max(length(alpha), length(pvalue))
  alpha   <- rep_len(alpha, n_combo)
  pvalue  <- rep_len(pvalue, n_combo)

  sd_obj <- .hdx_shared(df, state1, state2, times, replicates, alpha[1], pvalue[1], shared_data, "hdx_volcano")
  base_tbl <- sd_obj$data(withSelection = FALSE)

  ## positions (delta, -log10(p_value)) are the same for every combo - only
  ## is_sig (hence color) and the pvalue threshold line move - so each combo
  ## only needs its is_sig column, aligned to base_tbl's row order (the order
  ## actually plotted, which may come from an externally-supplied shared_data).
  raw_tables <- lapply(seq_len(n_combo), function(k) {
    .hdx_comparison_table(df, state1, state2, times, replicates, alpha[k], pvalue[k])
  })
  tables <- lapply(raw_tables, function(t) t[match(base_tbl$peptide_id, t$peptide_id), ])
  ## a peptide with identical replicate values in both states has an
  ## undefined Welch p-value (NA), so plotly's own build step silently drops
  ## that row from the trace (x/y incomplete) - match that filtering here too,
  ## or the marker.color vector below ends up one element longer than the
  ## rendered trace and every subsequent color shifts by one point.
  keep <- stats::complete.cases(base_tbl$delta, base_tbl$p_value)
  color_for <- function(tbl) ifelse(tbl$is_sig[keep], "firebrick", "grey60")
  vline_shapes <- function(pv) list(
    list(type = "line", x0 = 0, x1 = 0, xref = "x", y0 = -1, y1 = 1, yref = "paper",
         line = list(color = "grey80", dash = "dot")),
    list(type = "line", x0 = -1, x1 = 1, xref = "paper",
         y0 = -log10(pv), y1 = -log10(pv), yref = "y",
         line = list(color = "grey60", dash = "dash"))
  )

  p <- plotly::plot_ly(sd_obj, x = ~delta, y = ~-log10(p_value), text = ~hover, hoverinfo = "text",
                        type = "scatter", mode = "markers",
                        marker = list(size = 9, opacity = 0.75, color = color_for(tables[[1]]))) %>%
    plotly::layout(
      xaxis = list(title = "delta (state2 - state1)"),
      yaxis = list(title = "-log10(p-value)"),
      shapes = vline_shapes(pvalue[1])
    ) %>%
    plotly::highlight(on = "plotly_hover")

  if (n_combo > 1) {
    combo_lbl <- sprintf("alpha=%s, p=%s", alpha, pvalue)
    buttons <- lapply(seq_len(n_combo), function(k) {
      list(method = "update", label = combo_lbl[k],
           args = list(list(marker.color = list(color_for(tables[[k]]))),
                       list(shapes = vline_shapes(pvalue[k]))))
    })
    p <- plotly::layout(p, updatemenus = list(list(
      type = "dropdown", active = 0, x = 1.02, y = 1, xanchor = "left",
      pad = list(r = 10, t = 10), buttons = buttons
    )))
  }
  p
}

#' Interactive fissure_vent plot, linked to the volcano/Woods/table views
#'
#' Plotly version of [plot_fissure_vent()]: state1 vs state2 deuteration
#' plotted directly (no delta column), diagonal = no change, grey band =
#' the pooled Hageman & Weis critical interval, and each peptide drawn as
#' an ellipse built from that peptide's own Welch test - touching the
#' diagonal itself means that peptide's own p-value did not clear
#' `pvalue`. See `vignette("fissure_vent")` for the full explanation.
#'
#' @inheritParams interactive_volcano
#' @return a `plotly` htmlwidget.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' interactive_fissure_vent(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' }
#' @export
interactive_fissure_vent <- function(df, state1, state2, times, replicates = 3,
                                      alpha = 0.01, pvalue = 0.01, shared_data = NULL) {
  .hdx_require_interactive()
  n_combo <- max(length(alpha), length(pvalue))
  alpha   <- rep_len(alpha, n_combo)
  pvalue  <- rep_len(pvalue, n_combo)

  sd_obj <- .hdx_shared(df, state1, state2, times, replicates, alpha[1], pvalue[1], shared_data, "hdx_fissure")
  base_tbl <- sd_obj$data(withSelection = FALSE)

  ## positions (state1_val, state2_val) are the same for every combo - only
  ## is_sig/rx/ry (hence colors and ellipse shapes) and the band offset `ci`
  ## move - so each combo needs its own table (for rx/ry/is_sig) plus its own
  ## `ci`, aligned to base_tbl's row order (the order actually plotted, which
  ## may come from an externally-supplied shared_data).
  raw_tables <- lapply(seq_len(n_combo), function(k) {
    .hdx_comparison_table(df, state1, state2, times, replicates, alpha[k], pvalue[k])
  })
  cis <- vapply(raw_tables, function(t) attr(t, "ci"), numeric(1))
  tables <- lapply(raw_tables, function(t) t[match(base_tbl$peptide_id, t$peptide_id), ])

  rng <- range(c(base_tbl$state1_val, base_tbl$state2_val), na.rm = TRUE)
  pad <- diff(rng) * 0.08
  lims <- c(rng[1] - pad, rng[2] + pad)

  ## one ellipse per peptide, as an SVG path shape (48-gon approximation,
  ## axis-aligned - same construction as DrawEllipse() in the static plot)
  shapes_for <- function(tbl, ci) {
    ellipse_shapes <- lapply(seq_len(nrow(tbl)), function(i) {
      if (anyNA(c(tbl$state1_val[i], tbl$state2_val[i], tbl$rx[i], tbl$ry[i]))) return(NULL)
      theta <- seq(0, 2 * pi, length.out = 48)
      xs <- tbl$state1_val[i] + tbl$rx[i] * cos(theta)
      ys <- tbl$state2_val[i] + tbl$ry[i] * sin(theta)
      path <- paste0("M ", xs[1], ",", ys[1], " ",
                      paste(sprintf("L%s,%s", xs[-1], ys[-1]), collapse = " "), " Z")
      list(type = "path", path = path, xref = "x", yref = "y",
           fillcolor = if (isTRUE(tbl$is_sig[i])) "rgba(178,34,34,0.25)" else "rgba(120,120,120,0.15)",
           line = list(color = if (isTRUE(tbl$is_sig[i])) "firebrick" else "grey60", width = 1))
    })
    ellipse_shapes <- Filter(Negate(is.null), ellipse_shapes)
    band_shapes <- list(
      list(type = "line", x0 = lims[1], x1 = lims[2], y0 = lims[1], y1 = lims[2],
           line = list(color = "black")),
      list(type = "line", x0 = lims[1], x1 = lims[2], y0 = lims[1] + ci, y1 = lims[2] + ci,
           line = list(color = "grey70", dash = "dash")),
      list(type = "line", x0 = lims[1], x1 = lims[2], y0 = lims[1] - ci, y1 = lims[2] - ci,
           line = list(color = "grey70", dash = "dash"))
    )
    c(band_shapes, ellipse_shapes)
  }
  color_for <- function(tbl) ifelse(tbl$is_sig, "firebrick", "grey40")

  p <- plotly::plot_ly(sd_obj, x = ~state1_val, y = ~state2_val, text = ~hover, hoverinfo = "text",
                        type = "scatter", mode = "markers",
                        marker = list(size = 6, opacity = 0.9, color = color_for(tables[[1]]))) %>%
    plotly::layout(
      xaxis = list(title = paste("Deuteration,", attr(raw_tables[[1]], "state1")), range = lims),
      yaxis = list(title = paste("Deuteration,", attr(raw_tables[[1]], "state2")), range = lims),
      shapes = shapes_for(tables[[1]], cis[1])
    ) %>%
    plotly::highlight(on = "plotly_hover")

  if (n_combo > 1) {
    combo_lbl <- sprintf("alpha=%s, p=%s", alpha, pvalue)
    buttons <- lapply(seq_len(n_combo), function(k) {
      list(method = "update", label = combo_lbl[k],
           args = list(list(marker.color = list(color_for(tables[[k]]))),
                       list(shapes = shapes_for(tables[[k]], cis[k]))))
    })
    p <- plotly::layout(p, updatemenus = list(list(
      type = "dropdown", active = 0, x = 1.02, y = 1, xanchor = "left",
      pad = list(r = 10, t = 10), buttons = buttons
    )))
  }
  p
}

#' Interactive Woods plot, linked to the volcano/fissure_vent/table views
#'
#' Each peptide drawn as two horizontal segments across its `Start`-`End`
#' range, at height = state1's and state2's average deuteration
#' respectively (same layout convention as [woods_CI_plot()]), colored by
#' the same hybrid significance criterion as the other interactive plots.
#'
#' @inheritParams interactive_volcano
#' @return a `plotly` htmlwidget.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' interactive_woods(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' }
#' @export
interactive_woods <- function(df, state1, state2, times, replicates = 3,
                               alpha = 0.01, pvalue = 0.01, shared_data = NULL) {
  .hdx_require_interactive()
  sd_obj <- .hdx_shared(df, state1, state2, times, replicates, alpha, pvalue, shared_data, "hdx_woods")
  tbl <- sd_obj$data(withSelection = FALSE)

  seg <- function(val_col, label, col_sig, col_nonsig) {
    do.call(rbind, lapply(seq_len(nrow(tbl)), function(i) {
      data.frame(peptide_id = tbl$peptide_id[i], x = c(tbl$Start[i], tbl$End[i]),
                 y = tbl[[val_col]][i], state = label,
                 col = if (isTRUE(tbl$is_sig[i])) col_sig else col_nonsig,
                 hover = tbl$hover[i], stringsAsFactors = FALSE)
    }))
  }
  segs1 <- seg("state1_val", attr(tbl, "state1"), "firebrick", "lightpink")
  segs2 <- seg("state2_val", attr(tbl, "state2"), "darkblue", "lightblue")

  ## each state's segments become their own crosstalk::SharedData, sharing
  ## sd_obj's group name, so hovering here (or on the volcano/fissure_vent
  ## plots, or in the table) highlights the same peptide_id everywhere.
  ## Previously these were plain add_trace() calls with no crosstalk key at
  ## all, so the Woods plot never participated in the shared highlighting -
  ## it neither lit up when a peptide was hovered elsewhere, nor could it
  ## drive highlighting on the other panels itself.
  segs1_sd <- crosstalk::SharedData$new(segs1, key = ~peptide_id, group = sd_obj$groupName())
  segs2_sd <- crosstalk::SharedData$new(segs2, key = ~peptide_id, group = sd_obj$groupName())

  plotly::plot_ly() %>%
    plotly::add_trace(data = segs1_sd, x = ~x, y = ~y, split = ~peptide_id,
                       type = "scatter", mode = "lines", color = ~I(col),
                       line = list(width = 3), text = ~hover, hoverinfo = "text",
                       showlegend = FALSE) %>%
    plotly::add_trace(data = segs2_sd, x = ~x, y = ~y, split = ~peptide_id,
                       type = "scatter", mode = "lines", color = ~I(col),
                       line = list(width = 3), text = ~hover, hoverinfo = "text",
                       showlegend = FALSE) %>%
    plotly::layout(xaxis = list(title = "Residue"),
                   yaxis = list(title = "Deuteration"),
                   title = sprintf("%s (red) vs %s (blue)", attr(tbl, "state1"), attr(tbl, "state2"))) %>%
    plotly::highlight(on = "plotly_hover")
}

#' Interactive table of peptides, linked to the volcano/fissure_vent/Woods plots
#'
#' `DT::datatable()` of every peptide in the comparison (sortable/searchable),
#' with a `significant` column from the same hybrid criterion as the other
#' interactive plots. Selecting rows highlights the corresponding peptides
#' on any other plot built from the same `shared_data`.
#'
#' @inheritParams interactive_volcano
#' @return a `DT::datatable` htmlwidget.
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)
#'     && requireNamespace("DT", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' significant_peptides_table(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' }
#' @export
significant_peptides_table <- function(df, state1, state2, times, replicates = 3,
                                        alpha = 0.01, pvalue = 0.01, shared_data = NULL) {
  .hdx_require_interactive(need_dt = TRUE)
  sd_obj <- .hdx_shared(df, state1, state2, times, replicates, alpha, pvalue, shared_data, "hdx_table")

  ## build the display table and sort it as ONE data.frame, so peptide_id
  ## (the crosstalk join key) always travels with the row it identifies -
  ## sorting the key and the display columns separately previously left
  ## peptide_id pointing at the wrong row's Sequence/Start/End/Delta/p_value
  full <- sd_obj$data(withSelection = FALSE)
  full <- full[order(full$p_value), ]
  view <- data.frame(
    peptide_id = full$peptide_id, Sequence = full$Sequence,
    Start = full$Start, End = full$End,
    Delta = round(full$delta, 3), p_value = signif(full$p_value, 3),
    Significant = full$is_sig,
    stringsAsFactors = FALSE
  )

  sd_view <- crosstalk::SharedData$new(view, key = ~peptide_id, group = sd_obj$groupName())
  DT::datatable(sd_view, selection = "multiple", filter = "top",
                options = list(pageLength = 10, columnDefs = list(list(visible = FALSE, targets = 0))),
                rownames = FALSE)
}

#' Combined interactive explorer: volcano + fissure_vent + Woods + table
#'
#' Builds one shared dataset and lays out [interactive_volcano()],
#' [interactive_fissure_vent()], [interactive_woods()], and
#' [significant_peptides_table()] together, all linked via `crosstalk` -
#' selecting peptides in the table (or hovering over a point on any scatter
#' plot or Woods segment) highlights the same peptides everywhere. Requires
#' the `plotly`, `crosstalk`, and `DT` packages.
#'
#' @inheritParams interactive_volcano
#' @param alpha alpha for the pooled critical interval. Default `0.01`. May
#'   be a vector (paired element-wise with `pvalue`) - if so, the volcano
#'   and fissure_vent panels each get their own dropdown (see
#'   [interactive_volcano()]) for switching between the pre-computed
#'   `alpha`/`pvalue` combinations. This is a native per-plot `plotly`
#'   control, not a shared one: switching one panel's dropdown does not
#'   move the other panel's dropdown, the Woods plot, or the table - a
#'   single `htmlwidget`'s dropdown can't reach into a separate `htmlwidget`
#'   without custom JavaScript, which this package does not add. The Woods
#'   plot and table always reflect `alpha[1]`/`pvalue[1]` (the same
#'   combination the cross-widget hover-highlighting key is built from).
#' @param pvalue per-peptide significance cutoff. Default `0.01`. See `alpha`.
#' @return an `htmltools` tag list (a `crosstalk::bscols()` layout of
#'   several widgets, not a single htmlwidget) - viewable in the RStudio
#'   viewer, embeddable in an R Markdown document, or saved standalone with
#'   `htmltools::save_html()` (NOT `htmlwidgets::saveWidget()`, which only
#'   accepts a single widget).
#' @examples
#' if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)
#'     && requireNamespace("DT", quietly = TRUE)) {
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm)
#' interactive_fissure_explorer(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' interactive_fissure_explorer(a, state1 = "Unbound", state2 = "bound", times = "1800.00s",
#'                               alpha = c(0.01, 0.05), pvalue = c(0.01, 0.05))
#' }
#' @export
interactive_fissure_explorer <- function(df, state1, state2, times, replicates = 3,
                                          alpha = 0.01, pvalue = 0.01) {
  .hdx_require_interactive(need_dt = TRUE, need_htmltools = TRUE)
  sd_obj <- .hdx_shared(df, state1, state2, times, replicates, alpha[1], pvalue[1], NULL, "hdx_explorer")

  volcano <- interactive_volcano(df, state1, state2, times, replicates, alpha, pvalue, shared_data = sd_obj)
  fissure <- interactive_fissure_vent(df, state1, state2, times, replicates, alpha, pvalue, shared_data = sd_obj)
  woods <- interactive_woods(df, state1, state2, times, replicates, alpha[1], pvalue[1], shared_data = sd_obj)
  tbl <- significant_peptides_table(df, state1, state2, times, replicates, alpha[1], pvalue[1], shared_data = sd_obj)

  htmltools::browsable(htmltools::tagList(
    crosstalk::bscols(widths = c(6, 6), volcano, fissure),
    crosstalk::bscols(widths = 12, woods),
    crosstalk::bscols(widths = 12, tbl)
  ))
}
