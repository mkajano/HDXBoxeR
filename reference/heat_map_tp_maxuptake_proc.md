# Preparatory function for heat map of maximum procent deuteration per residue.

Returns heat map

## Usage

``` r
heat_map_tp_maxuptake_proc(
  df,
  dfup,
  pv,
  sd,
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  replicates = 3,
  pv_cutoff
)
```

## Arguments

- df:

  average data frame for procent deuteration. Generated using
  ave_timepoint() function.

- dfup:

  average data frame for deuteration uptake. Generated using
  ave_timepoint() function.

- pv:

  pvalues dataframes calculated using pv_timepoint() function

- sd:

  standard deviation data.frame generated using sd_timepoint function

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- replicates:

  number of replicates in sample. Default set to 3.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

Maximum uptake heat map for timepoints
