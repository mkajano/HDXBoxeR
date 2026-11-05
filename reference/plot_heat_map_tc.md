# Plots heat maps for time courses.

Returns heat map on timecourses with raw data.

## Usage

``` r
plot_heat_map_tc(
  df,
  replicates = 3,
  mar_x = 3.5,
  ranges = c(-Inf, seq(0, 100, by = 10), Inf)
)
```

## Arguments

- df:

  output from function output_tc

- replicates:

  number of replicates in sample. Default set to 3.

- mar_x:

  margin x width. Default=3.5

- ranges:

  ranges for coloring scheme. Default set to c(seq(0, 100, by=10), Inf)

## Value

heat map for time courses

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tc(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
plot_heat_map_tc(df=a, replicates=3, ranges=c(seq(0, 100, by=5), Inf))
#> values range from  0 to 10.89 %
#> t3.00s range 0 8.54
#> t60.00s range 0 9.87
#> t1800.00s range 0 10.43
#> t72000.00s range 0 10.89
#> values range from  0 to 10.4 %
#> t3.00s range 0 7.73
#> t60.00s range 0 8.98
#> t1800.00s range 0 9.61
#> t72000.00s range 0 10.4

#> NULL
plot_heat_map_tc(df=a)
#> values range from  0 to 10.89 %
#> t3.00s range 0 8.54
#> t60.00s range 0 9.87
#> t1800.00s range 0 10.43
#> t72000.00s range 0 10.89
#> values range from  0 to 10.4 %
#> t3.00s range 0 7.73
#> t60.00s range 0 8.98
#> t1800.00s range 0 9.61
#> t72000.00s range 0 10.4

#> NULL
```
