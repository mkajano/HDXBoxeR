# Returns volcano plots for timepoints in the data frame

Returns volcano plots for each peptide. Critical interval is calculated
according to \#' Reliable Identification of Significant Differences in
Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a
Hybrid Significance Testing Approach Tyler S. Hageman and David D. Weis
Analytical Chemistry 2019 91 (13), 8008-8016 DOI:
10.1021/acs.analchem.9b01325 pvalues calculated using Welch t-test. Each
state is colored with its own critical interval (vs the first state).

## Usage

``` r
plots_vol_tp(
  df,
  replicates = 3,
  pvalue = 0.01,
  alpha = 0.01,
  col,
  p_adjust = "none",
  pv_cutoff,
  cola
)
```

## Arguments

- df:

  output from functions output_tp

- replicates:

  number of replicates in set as default set to 3.

- pvalue:

  p-value cutoff here set up to 0.01

- alpha:

  cutoff for Critical interval. Default=0.01

- col:

  color pallette for different Protein States. As default Paired
  pallette from color.Brewer is used.

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- pv_cutoff:

  deprecated, use `pvalue`.

- cola:

  deprecated, use `col`.

## Value

volcano plots

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
plots_vol_tp(df=a, replicates=3, col=c(1:4), pvalue=0.01, alpha=0.01 )
#> CI @ 3.00s : 0.35
#> CI @ 60.00s : 0.27
#> CI @ 1800.00s : 0.24
#> CI @ 72000.00s : 0.2

plots_vol_tp(df=a, pvalue=0.05)
#> CI @ 3.00s : 0.35
#> CI @ 60.00s : 0.27
#> CI @ 1800.00s : 0.24
#> CI @ 72000.00s : 0.2
```
