# Generates average deuteration plot for the time-course.

Returns plots with average deuteration at each peptide.

## Usage

``` r
plots_av_tcourse(df, replicates = 3, col, cola)
```

## Arguments

- df:

  output from functions output_tc (uptake or percent = TRUE).

- replicates:

  number of replicates in set as default set to 3.

- col:

  color pallette for different Protein States. As default Paired
  pallette from RColorBrewer is used.

- cola:

  deprecated, use `col`.

## Value

average deuteration plots

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tc(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
plots_av_tcourse(df=a, replicates=3, col=c(1:4))

plots_av_tcourse(df=a)
```
