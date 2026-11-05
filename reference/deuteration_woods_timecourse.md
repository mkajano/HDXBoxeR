# Return woods plots for the timecourse

All the peptides are plotted based on their uptake.

## Usage

``` r
deuteration_woods_timecourse(
  input_data,
  states,
  replicates = 3,
  ylim = c(0, 120),
  ...
)
```

## Arguments

- input_data:

  output from function output_tc(..., percent=TRUE)

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

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tc(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
deuteration_woods_timecourse(a)

# }
```
