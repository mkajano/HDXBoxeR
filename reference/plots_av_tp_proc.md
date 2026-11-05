# Returns average procent deuteration plot for time points

Returns plots with average procent deuteration at each peptide.

## Usage

``` r
plots_av_tp_proc(df, replicates = 3, col, cola)
```

## Arguments

- df:

  output from functions output_tp_proc.

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
a<- output_tp(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
plots_av_tp_proc(df=a, replicates=3, col=c(1:4))

plots_av_tp_proc(df=a)
```
