# Legend for the heatmaps.Extracts names from data.frame

Returns names for legend for the heatmaps

## Usage

``` r
legend_heat_map_tp(df)
```

## Arguments

- df:

  average data frame. Generated using ave_timepoint() function.

## Value

legend for the heatmap

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
legend_heat_map_tp(df=a)
#> Error in axis(2, at = 0.5:(dim(df)[2] - 7.5), labels = nm1, las = 2, line = ,     cex.axis = 0.7): plot.new has not been called yet
```
