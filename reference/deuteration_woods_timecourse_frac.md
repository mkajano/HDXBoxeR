# Return woods plots for the timecourse

All the peptides are plotted based on their uptake.

## Usage

``` r
deuteration_woods_timecourse_frac(
  input_data,
  states,
  replicates = 3,
  ylim = c(0, 1.2),
  ...
)
```

## Arguments

- input_data:

  output from function output_tc() with uptake in Da (it is divided by
  the number of exchangeable amides, see nb_exch_deut())

- states:

  states, if missing all states used

- replicates:

  replicates

- ylim:

  y axis limits

- ...:

  other parameters

## Value

Woods plots for the timecourse
