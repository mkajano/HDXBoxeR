# Legend for [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)'s color scale

Thin wrapper around
[`pallette_ll()`](https://mkajano.github.io/HDXBoxeR/reference/pallette_ll.md) -
draws the same red/blue bins used by
[`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)
as a base-R legend, since the color scale isn't otherwise labeled inside
the 3D `NGLVieweR` widget itself.

## Usage

``` r
plot_hdx_structure_legend(
  ranges,
  col = "grey80",
  unit_label = "",
  colors_initial
)
```

## Arguments

- ranges:

  the same `ranges` vector passed to
  [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)
  (or, if that call used the default `NULL`, the auto-scaled vector -
  see
  [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)'s
  `@details` on reconstructing it).

- col:

  the same `col` passed to
  [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md).
  Default `"grey80"`.

- unit_label:

  axis label for the bins, e.g. `"Da"` or `"%"`. Default `""`.

- colors_initial:

  deprecated, use `col`.

## Value

draws a base-R legend plot; returns `NULL` invisibly.

## See also

[`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)

## Examples

``` r
plot_hdx_structure_legend(ranges = c(-Inf, -20, -10, 0, 10, 20, Inf), unit_label = "Da")
```
