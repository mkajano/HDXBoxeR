# Plots heat maps for maximum uptake per residue.

Returns heat map with maximum uptake per residue.

## Usage

``` r
plot_heat_map_max_uptake_tp(
  df,
  replicates = 3,
  mar_x = 3.5,
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  pv_cutoff
)
```

## Arguments

- df:

  average data frame. Generated using ave_timepoint() function.

- replicates:

  number of replicates in sample. Default set to 3.

- mar_x:

  margin x width. Default=3.5

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

heat map for maximum uptake per residue

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
plot_heat_map_max_uptake_tp(df=a, replicates=3, alpha=0.01,pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )
#> For timepoint 3.00s
#> values range from  0 to 5.68 %
#> Unbound bound range 0 5.68
#> For timepoint 60.00s
#> values range from  0 to 4.04 %
#> Unbound bound range 0 4.04
#> For timepoint 1800.00s
#> values range from  0 to 1.96 %
#> Unbound bound range 0 1.96
#> For timepoint 72000.00s
#> values range from  -0.99 to 6.01 %
#> Unbound bound range -0.99 6.01

#> NULL
plot_heat_map_max_uptake_tp(df=a)
#> For timepoint 3.00s
#> values range from  0 to 5.68 %
#> Unbound bound range 0 5.68
#> For timepoint 60.00s
#> values range from  0 to 4.04 %
#> Unbound bound range 0 4.04
#> For timepoint 1800.00s
#> values range from  0 to 1.96 %
#> Unbound bound range 0 1.96
#> For timepoint 72000.00s
#> values range from  -0.99 to 6.01 %
#> Unbound bound range -0.99 6.01

#> NULL
```
