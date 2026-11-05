# Prepares function for Critical interval for timecourses

Preparatory function

## Usage

``` r
prep_timecourse_plot_sd(
  control_df_up,
  variant_df_up,
  replicates = 3,
  alpha = 0.01
)
```

## Arguments

- control_df_up:

  dataframe of control

- variant_df_up:

  dataframe for variant

- replicates:

  number of replicates. Default set to 3.

- alpha:

  cut off of pvalue used in calculation of critical interval. Default
  set to 0.01

## Value

Critial interval for all sets
