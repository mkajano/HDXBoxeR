# Return woods plots for the timepoints

All the peptides are plotted based on their uptake.

## Usage

``` r
deuteration_woods_timepoints(
  input_data,
  times,
  replicates = 3,
  col = NA,
  ylim = c(0, 120),
  ...,
  cola
)
```

## Arguments

- input_data:

  output from function output_tp(..., percent=TRUE)

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

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
deuteration_woods_timepoints(a[1:12,])

# }
```
