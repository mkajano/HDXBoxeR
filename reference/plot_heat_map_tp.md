# Plots heat maps for significant peptides.

Returns heat map with average values for significant uptake per residue.

## Usage

``` r
plot_heat_map_tp(
  df,
  mar_x = 3.5,
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  replicates = 3,
  p_adjust = "none",
  pv_cutoff
)
```

## Arguments

- df:

  average data frame. Generated using ave_timepoint() function.

- mar_x:

  margin x width. Default=3.5

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- replicates:

  number of replicates in sample. Default set to 3.

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

heat map for average uptake per residue for significant peptides.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
plot_heat_map_tp(df=a, replicates=3, alpha=0.01,pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )
#> For time point 3.00s
#> values range from  0 to 10.35 %
#> Unbound bound range 0 10.35
#> For time point 60.00s
#> values range from  0 to 6.71 %
#> Unbound bound range 0 6.71
#> For time point 1800.00s
#> values range from  0 to 3.11 %
#> Unbound bound range 0 3.11
#> For time point 72000.00s
#> values range from  0 to 2.81 %
#> Unbound bound range 0 2.81

#> NULL
plot_heat_map_tp(df=a)
#> For time point 3.00s
#> values range from  0 to 10.35 %
#> Unbound bound range 0 10.35
#> For time point 60.00s
#> values range from  0 to 6.71 %
#> Unbound bound range 0 6.71
#> For time point 1800.00s
#> values range from  0 to 3.11 %
#> Unbound bound range 0 3.11
#> For time point 72000.00s
#> values range from  0 to 2.81 %
#> Unbound bound range 0 2.81

#> NULL
```
