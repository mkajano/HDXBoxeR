# Plots boxplots for all the averages in the set

Returns boxplots to compare sets between each other

## Usage

``` r
boxplot_tp(df, replicates = 3, ...)
```

## Arguments

- df:

  average data frame. Generated using ave_timepoint() function.

- replicates:

  number of replicates in sample. Default set to 3.

- ...:

  inherited boxplot parameters

## Value

boxplots for average deuterium uptake per set.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
boxplot_tp(df=a, replicates=3)
```
