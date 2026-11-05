# Uptake plots

Uptake plots per peptide

## Usage

``` r
uptake_plots(
  input_data,
  timepoints,
  replicates = 3,
  col = NA,
  seq_match = TRUE,
  cola
)
```

## Arguments

- input_data:

  output from function output_tc(..., percent=TRUE)

- timepoints:

  the labeling times

- replicates:

  replicates

- col:

  colors, default NA

- seq_match:

  Flag TRUE or FALSE, default TRUE, match sequence of the protein states

- cola:

  deprecated, use `col`.

## Value

Uptake plots

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
x=c(3,60, 1800, 72000)
uptake_plots(a, x)







































# }
```
