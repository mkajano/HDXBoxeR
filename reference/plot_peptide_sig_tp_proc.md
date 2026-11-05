# Draws peptides with significant difefrences between sets.

Returns plot where significant peptides are colored in blue-red scheme.

## Usage

``` r
plot_peptide_sig_tp_proc(
  input_proc,
  input_up,
  nb_pep_row = 100,
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

- nb_pep_row:

  number of peptides in each row. Plotting parameter. Default set to
  100.

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

plot with peptides which are significantly different between sets.

## Details

Percent-deuteration counterpart of
[`plot_peptide_sig_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_sig_tp.md):
significance is tested on the uptake table `input_up` (Da, the
package-wide rule), while peptides are colored by their difference in
percent deuteration, first state minus each state (from `input_proc`);
[`plot_peptide_sig_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_sig_tp.md)
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
plot_peptide_sig_tp_proc(input_proc=a_proc, input_up=a_up, replicates=3,
alpha=0.01, pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf), nb_pep_row=40 )
#> For time point 3.00s
#> values range from  -0.99 to 3.56 %
#> Unbound bound range -0.99 3.56

#> For time point 60.00s
#> values range from  -1.33 to 2.54 %
#> Unbound bound range -1.33 2.54

#> For time point 1800.00s
#> values range from  -1.12 to 2.2 %
#> Unbound bound range -1.12 2.2

#> For time point 72000.00s
#> values range from  -2.98 to 0.7 %
#> Unbound bound range -2.98 0.7
```
