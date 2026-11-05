# Return woods plots for the timepoints

All the peptides are plotted based on their uptake.

## Usage

``` r
deuteration_woods_timepoints_frac(
  input_data,
  times,
  replicates = 3,
  col = NA,
  ylim = c(0, 1.1),
  ...,
  cola
)
```

## Arguments

- input_data:

  output from function output_tp() with uptake in Da (it is divided by
  the number of exchangeable amides, see nb_exch_deut())

- times:

  Deuteration times, if missing all deuteration times used

- replicates:

  replicates

- col:

  colors, default NA

- ylim:

  y axis limits

- ...:

  other parameters

- cola:

  deprecated, use `col`.

## Value

Woods plots for the timepoints
