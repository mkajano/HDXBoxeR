# pvalue calculation between two sets of the data at certain timepoint

Preparatory function for calculation of pvalue between sets.

## Usage

``` r
pv_timecourse(df_c, df_v, replicates = 3, p_adjust = "none")
```

## Arguments

- df_c:

  dataframe of control

- df_v:

  dataframe for variant

- replicates:

  number of replicates. Default set to 3.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default) or e.g. `"BH"`, applied across the peptides of the
  comparison, separately for each time point.

## Value

pvalue comparisons between two sets.
