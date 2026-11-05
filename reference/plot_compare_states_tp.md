# Plot per-peptide deltas from `compare_states_tp()`, by state pair and region

[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
itself has no plotting companion - unlike
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
(paired with
[`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md))
its peptide-by-pair result table was only ever inspected as a
data.frame. This draws one panel per region (or a single "all" panel if
`regions` wasn't used), with every matched peptide's `delta` shown as a
box, violin, or plain jittered-point strip - one position per state
pair - so a whole domain's behavior for a given comparison is visible at
a glance, for any number of states.

## Usage

``` r
plot_compare_states_tp(
  res,
  plot_type = c("jitter", "box", "violin"),
  sig_col = "firebrick",
  nonsig_col = "grey55",
  summary_col = "steelblue",
  tint = 0.6,
  point_cex = 0.9,
  jitter_width = 0.15,
  show_zero = TRUE,
  bg = "grey95",
  grid_col = "white",
  show_grid = TRUE,
  cex_annot = 0.9,
  font_family = ""
)
```

## Arguments

- res:

  the data.frame returned by
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md).

- plot_type:

  one of `"jitter"` (every peptide as a point, no summary shape),
  `"box"`, or `"violin"`. Default `"jitter"`.

- sig_col:

  color for peptides with `is_sig == TRUE`. Default `"firebrick"`.

- nonsig_col:

  color for peptides with `is_sig == FALSE`. Default `"grey55"`.

- summary_col:

  color for the box/violin summary shape (drawn once per pair, pooling
  significant and non-significant peptides together). Default
  `"steelblue"`.

- tint:

  how much lighter than `summary_col` the box/violin fill is, `0` (no
  lightening) to `1` (white). Default `0.6`.

- point_cex:

  size of the jittered points. Default `0.9`.

- jitter_width:

  horizontal jitter width for points, as a fraction of pair spacing.
  Default `0.15`.

- show_zero:

  logical; draw a dashed horizontal reference line at `delta = 0`.
  Default `TRUE`.

- bg:

  panel background color. Default `"grey95"`; use `"white"` for a plain
  white panel.

- grid_col:

  gridline color drawn over `bg`. Default `"white"`.

- show_grid:

  logical; draw the shaded background/gridlines at all. Default `TRUE`.

- cex_annot:

  text size for the "n significant / n tested" annotation under each
  region's title. Default `0.9`.

- font_family:

  font family passed to text-drawing calls. Default `""` (device
  default).

## Value

draws a one-panel-per-region base-R plot; returns `NULL` invisibly.

## Showing every peptide, not only the significant ones

Every peptide tested is drawn, regardless of `is_sig` - significant ones
are simply colored differently (`sig_col` vs `nonsig_col`) rather than
filtered out. Looking only at a list of "significant peptides" hides
exactly the pattern
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)'s
`frac_peptides_significant` diagnostic exists to catch: a region can
contain a real, consistent shift in most of its peptides while still
having only a few (or zero) individually clear the strict Hageman & Weis
gate, or vice versa - one dramatic outlier peptide surrounded by
otherwise flat ones. Seeing the whole distribution, not just the subset
that cleared a cutoff, is the point of this plot.

## On color and visibility

Box/violin fills are opaque colors (lightened by blending toward white
in RGB space via `tint`, not by alpha transparency) deliberately - some
graphics devices render semi-transparent fills as invisible or fully
solid instead of blending them, so visibility here never depends on
device compositing.

## See also

[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md),
[`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md)

## Examples

``` r
# \donttest{
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
# the shipped example only has two states; compare_states_tp() needs 3+.
# See vignette("multistate") for a worked example with a synthetic third state.
# }
```
