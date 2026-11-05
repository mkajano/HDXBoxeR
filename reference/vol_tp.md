# Preparatory function for volcano plot

Returns volcano plots

## Usage

``` r
vol_tp(df1, pv, CI, pvalue = 0.01, col, pv_cutoff, cola)
```

## Arguments

- df1:

  differences in averages data.frame calculated using diff_ave function

- pv:

  pvalues dataframes calculated using pv_timepoint function

- CI:

  critical interval: one value per compared state (as from
  [`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md)),
  or a single value used for all.

- pvalue:

  p-value cutoff here set up to 0.01

- col:

  color pallette for different Protein States. As default Paired
  pallette from color.Brewer is used.

- pv_cutoff:

  deprecated, use `pvalue`.

- cola:

  deprecated, use `col`.

## Value

volcano plots
