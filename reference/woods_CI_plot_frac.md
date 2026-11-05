# Returns a woods plot for comparisons of the timepoints samples

Modification of butterfly plot. x axis residues. y axis % deuteration
for Peptides are compared between the sets for the significance change
between sets. If there is significant change beteween sets peptides are
plotted for all timepoints. Significanty different timepoints for the
peptides are colored. Peptides ranges are plotted as a line at
corresponding % deuteration values.

## Usage

``` r
woods_CI_plot_frac(
  th,
  replicates = 3,
  pvalue = 0.01,
  states,
  alpha = 0.01,
  ylim = c(0, 1.2),
  p_adjust = "none",
  ...,
  pv_cutoff
)
```

## Arguments

- th:

  output of output_tc(). Raw data for uptake deuteration for time
  courses

- replicates:

  number of replicates in sample. Default set to 3.

- pvalue:

  p-value cutoff here set up to 0.01

- states:

  Protein states from the set. As default all states are chosen.

- alpha:

  critical interval, to have more restrictive use lower values,
  default=0.01

- ylim:

  y axis limit

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- ...:

  other variables

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

Woods plots with chosen statistically different peptides

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tc(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
b<-output_tc(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
woods_CI_plot_frac(th=a, pvalue = 0.01, alpha = 0.01, replicates=3)

# }
```
