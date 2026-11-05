# Significant peptide plots.

Returns plot where significant peptides are colored in blue-red scheme.

## Usage

``` r
plot_peptide_sig_tp(
  df1,
  replicates = 3,
  nb_pep_row = 100,
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  pv_cutoff
)
```

## Arguments

- df1:

  average data frame. Generated using ave_timepoint() function.

- replicates:

  number of replicates in sample. Default set to 3.

- nb_pep_row:

  number of peptides in each row. Plotting parameter. Default set to
  100.

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

plot with peptides which are significantly different between sets.
