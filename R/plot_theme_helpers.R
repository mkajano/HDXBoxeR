## internal plotting-style helpers - not exported. Shared by
## plot_peptide_coverage_map(), plot_coverage_diagnostics(), and
## plot_compare_states_unmatched() so the three "ANOVA workflow" plots read
## as one visual system instead of three independently-styled base R plots.
## Deliberately base graphics only (no ggplot2 dependency) - a shaded panel
## background plus white gridlines mimics the "theme_bw"/"theme_minimal"
## look inside plot.window() coordinates, drawn AFTER plot(..., type = "n")
## and BEFORE the data layer, matching how ggplot2 itself layers panel,
## grid, then geoms.

## Resolves a `palette` argument to exactly n colors, consistently across
## every plot in this file group. `palette` is either the name of an
## RColorBrewer qualitative palette (checked against brewer.pal.info, not
## just "is it length 1" - a single explicit color name like "steelblue" is
## also length 1 and must NOT be sent to brewer.pal(), which would error)
## or an explicit vector of one or more colors, recycled if only one is
## given.
.hdx_palette <- function(n, palette = "Set2") {
  is_brewer_name <- length(palette) == 1 && palette %in% rownames(RColorBrewer::brewer.pal.info)
  if (!is_brewer_name) {
    if (length(palette) == 1) palette <- rep(palette, n)
    if (length(palette) < n) stop(sprintf("palette has %d colors but %d are needed", length(palette), n))
    return(palette[seq_len(n)])
  }
  pal_n <- max(3, n)
  RColorBrewer::brewer.pal(pal_n, palette)[seq_len(n)]
}

## Draws a shaded panel background and gridlines inside the CURRENT plot
## region (call right after plot(..., type = "n"), before adding data).
## `xlim`/`ylim` should be the same limits already passed to plot().
.hdx_panel <- function(xlim, ylim, bg = "grey95", grid_col = "white", grid_lwd = 1.2,
                        nx = NULL, ny = NULL) {
  graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2], col = bg, border = NA)
  if (!is.null(nx) && nx > 0) graphics::abline(v = pretty(xlim, nx), col = grid_col, lwd = grid_lwd)
  if (!is.null(ny) && ny > 0) graphics::abline(h = pretty(ylim, ny), col = grid_col, lwd = grid_lwd)
  if (is.null(nx)) graphics::abline(v = graphics::axTicks(1), col = grid_col, lwd = grid_lwd)
  if (is.null(ny)) graphics::abline(h = graphics::axTicks(2), col = grid_col, lwd = grid_lwd)
  graphics::box(col = "grey40")
}

## Lightens a color by blending it toward white in RGB space and returns an
## OPAQUE result - deliberately not grDevices::adjustcolor(alpha.f=...).
## True alpha transparency is compositing performed by the graphics DEVICE,
## and several base R devices (notably the Windows native `windows()`
## device RStudio's Plots pane can fall back to) render semi-transparent
## fills and points as fully invisible or fully solid instead of blending -
## a real, reported failure mode ("I cannot see the points, something is
## wrong with the opacity") even though the identical code renders
## correctly through a bitmap device like png()/cairo. Blending in RGB
## space ourselves produces a plain solid color, so visibility can never
## depend on which device is active.
## `amount` in [0, 1]: 0 = unchanged color, 1 = white.
.hdx_tint <- function(col, amount = 0.5) {
  rgb_col <- grDevices::col2rgb(col) / 255
  mixed <- rgb_col * (1 - amount) + amount
  grDevices::rgb(mixed[1, ], mixed[2, ], mixed[3, ])
}

## Draws one vertical violin (mirrored kernel density) centered at
## `x_center`, base graphics only (no dependency). Falls back to plotting
## the raw points if there are too few values to estimate a density.
.hdx_violin <- function(x_center, values, half_width = 0.38, col, border = "grey25", lwd = 1) {
  values <- values[!is.na(values)]
  if (length(values) < 3 || stats::sd(values) == 0) {
    graphics::points(rep(x_center, length(values)), values, pch = 16, col = col)
    return(invisible(NULL))
  }
  d <- stats::density(values, n = 128)
  w <- d$y / max(d$y) * half_width
  graphics::polygon(c(x_center - w, rev(x_center + w)), c(d$x, rev(d$x)),
                     col = col, border = border, lwd = lwd)
}

## Bold title + grey subtitle line, positioned consistently across plots.
.hdx_title <- function(title, subtitle = NULL, cex_title = 1.15, cex_sub = 0.9,
                        title_col = "grey15", sub_col = "grey40") {
  if (!is.null(subtitle)) {
    graphics::mtext(subtitle, side = 3, line = 0.3, adj = 0, cex = cex_sub, col = sub_col)
    graphics::mtext(title, side = 3, line = 1.5, adj = 0, cex = cex_title, font = 2, col = title_col)
  } else {
    graphics::mtext(title, side = 3, line = 0.5, adj = 0, cex = cex_title, font = 2, col = title_col)
  }
}
