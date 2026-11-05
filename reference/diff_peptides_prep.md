# Peptides significantly different between states

Returns list of the dataframes with significantly different peptides

## Usage

``` r
diff_peptides_prep(
  df,
  replicates = 3,
  pvalue = 0.01,
  alpha = 0.01,
  p_adjust = "none",
  pv_cutoff
)
```

## Arguments

- df:

  differences in averages data.frame calculated using diff_ave function

- replicates:

  number of replicates

- pvalue:

  p-value cutoff here set up to 0.01

- alpha:

  critical interval, here is multiple sets are using maximum CI is used.

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

list of the dataframes
