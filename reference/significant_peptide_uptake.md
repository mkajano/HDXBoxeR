# Function returns which peptides are significantly based of pv_cutoff and Critial interval

Returns data frame with significant peptides.

## Usage

``` r
significant_peptide_uptake(
  df_av,
  pv,
  sd,
  alpha = 0.01,
  pvalue = 0.01,
  replicates = 3,
  pv_cutoff
)
```

## Arguments

- df_av:

  data.frame with averages created using ave_timepoint() function

- pv:

  data.frame with pvalues created using pv_timepoint() function

- sd:

  data.frame with standard deviations created using sd_timepoint()
  function

- alpha:

  cutoff for Critical interval. Default=0.01

- pvalue:

  cuttoff for pval. Default=0.01

- replicates:

  number of replicates as default set to 3.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

ranges per set
