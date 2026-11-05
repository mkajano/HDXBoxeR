# Preparatory function for showing peptides with significant differences between sets.

Returns plot where significantly different peptides are colored in
blue-red scheme.

## Usage

``` r
peptide_pv_tp_proc(
  df,
  dfup,
  pv,
  sd,
  nb_row = 100,
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

- nb_row:

  number of peptides in each row. Plotting parameter.

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for Critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- replicates:

  number of replicates in sample. Default set to 3.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

plot with peptides which are significantly different between sets.
