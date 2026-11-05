# Calculation of pvalue between first protein state and any other state from all_states file

Welch t-test of every state (or time) against the first one in `df`, per
peptide. Missing replicates are ignored; the p-value is `NA` when either
side has fewer than two values or the test is undefined (e.g. all
replicates identical).

## Usage

``` r
pv_timepoint(df, replicates = 3, p_adjust = "none")
```

## Arguments

- df:

  output from functions output_tp or output_tc.

- replicates:

  only used for tables whose value columns are not named
  `<state>_..._<replicate>`: the number of replicates per block. Default
  3.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw Welch p-values) or e.g. `"BH"`. P-values are adjusted
  across the peptides of each compared state, separately for each
  labeling time (the same family as everywhere in the package).

## Value

Data.frame with p-values

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
pv<-pv_timepoint(df=a)
pv_bh<-pv_timepoint(df=a, p_adjust="BH")
```
