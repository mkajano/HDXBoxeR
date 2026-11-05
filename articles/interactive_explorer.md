# Interactive Volcano, Fissure, and Woods Plots

The interactive plots are for exploring the data. Hover over a peptide
to see its sequence, uptake, difference and p-value. Plots built
together are linked: a peptide selected in one is highlighted in all of
them. For figures in a paper, use the static plots
([`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
[`plots_vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_vol_tp.md)).

The plots need the plotly, crosstalk and DT packages. They are optional
and are not installed with HDXBoxeR.

| Function | Shows |
|----|----|
| [`interactive_volcano()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_volcano.md) | volcano plot |
| [`interactive_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_vent.md) | fissure vent plot ([`vignette("fissure_vent")`](https://mkajano.github.io/HDXBoxeR/articles/fissure_vent.md)) |
| [`interactive_woods()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_woods.md) | Woods plot |
| [`significant_peptides_table()`](https://mkajano.github.io/HDXBoxeR/reference/significant_peptides_table.md) | table of all peptides |
| [`interactive_fissure_explorer()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_explorer.md) | all four together, linked |

The significance test is the same as in the rest of the package.

## Example

``` r

library(HDXBoxeR)
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm)
```

One plot:

``` r

interactive_fissure_vent(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
```

Give several `alpha` and `pvalue` values to switch between them with a
dropdown (top right):

``` r

interactive_fissure_vent(a, state1 = "Unbound", state2 = "bound", times = "1800.00s",
                         alpha = c(0.01, 0.05), pvalue = c(0.01, 0.05))
```

`alpha` and `pvalue` are paired: the first with the first, the second
with the second.

All plots together. Select rows in the table, or hover over a peptide,
to highlight it everywhere:

``` r

interactive_fissure_explorer(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
```

## Saving

Save one plot with
[`htmlwidgets::saveWidget()`](https://rdrr.io/pkg/htmlwidgets/man/saveWidget.html).
The explorer is a layout of several plots, so save it with
[`htmltools::save_html()`](https://rstudio.github.io/htmltools/reference/save_html.html):

``` r

explorer <- interactive_fissure_explorer(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
htmltools::save_html(explorer, "fissure_explorer.html")
```

The saved file opens in any browser, without R.

## Good to know

- One time point per call.
- `p_adjust = "BH"` corrects the p-values for multiple testing.

## How this connects to the rest of the package

| Interactive | Static version | Explained in |
|----|----|----|
| [`interactive_volcano()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_volcano.md) | [`plots_vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_vol_tp.md) | [`vignette("HDXBoxeR")`](https://mkajano.github.io/HDXBoxeR/articles/HDXBoxeR.md) |
| [`interactive_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_vent.md) | [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md) | [`vignette("fissure_vent")`](https://mkajano.github.io/HDXBoxeR/articles/fissure_vent.md) |
| [`interactive_woods()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_woods.md) | [`woods_CI_plot()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot.md) | [`vignette("HDXBoxeR")`](https://mkajano.github.io/HDXBoxeR/articles/HDXBoxeR.md) |
| [`significant_peptides_table()`](https://mkajano.github.io/HDXBoxeR/reference/significant_peptides_table.md) | [`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md), [`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md) | [`vignette("HDXBoxeR")`](https://mkajano.github.io/HDXBoxeR/articles/HDXBoxeR.md) |

To show results on a structure, see
[`vignette("hdx_structure")`](https://mkajano.github.io/HDXBoxeR/articles/hdx_structure.md).
