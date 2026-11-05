# Plots heat maps for maximum procent deuteration per residue.

Returns heat map with maximum precent_deuteration per residue.

## Usage

``` r
plot_heat_map_max_uptake_tp_proc(
  input_proc,
  input_up,
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

- input_proc:

  Dataframe with organized procent deuteration data. Input generated
  using output_tp(percent = TRUE).

- input_up:

  Dataframe with organized deuteration uptake. Input generated using
  output_tp() function.

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

## Details

Percent-deuteration counterpart of
[`plot_heat_map_max_uptake_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_max_uptake_tp.md):
significance is tested on the uptake table `input_up` (Da, the
package-wide rule), while residues are colored by the largest
significant difference in percent deuteration, first state minus each
state (from `input_proc`) among the peptides covering them;
[`plot_heat_map_max_uptake_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_max_uptake_tp.md)
uses the relative change in uptake.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a_up<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
a_proc<- output_tp(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
plot_heat_map_max_uptake_tp_proc(input_proc=a_proc, input_up=a_up, replicates=3, pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )
#> For time point 3.00s
#> values range from  -9.86 to 35.61 %
#> Unbound bound range -9.86 35.61
#> For time point 60.00s
#> values range from  -13.3 to 25.36 %
#> Unbound bound range -13.3 25.36
#> For time point 1800.00s
#> values range from  -11.17 to 21.97 %
#> Unbound bound range -11.17 21.97
#> For time point 72000.00s
#> values range from  -29.81 to 7 %
#> Unbound bound range -29.81 7

#> NULL
plot_heat_map_max_uptake_tp_proc(input_proc=a_proc, input_up=a_up)
#> For time point 3.00s
#> values range from  -9.86 to 35.61 %
#> Unbound bound range -9.86 35.61
#> For time point 60.00s
#> values range from  -13.3 to 25.36 %
#> Unbound bound range -13.3 25.36
#> For time point 1800.00s
#> values range from  -11.17 to 21.97 %
#> Unbound bound range -11.17 21.97
#> For time point 72000.00s
#> values range from  -29.81 to 7 %
#> Unbound bound range -29.81 7

#> NULL
```
