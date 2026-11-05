# Returns indexes for peptides with significant difference between two sets

Function to help decide which peptides will be drawn on Robot plots.

## Usage

``` r
robot_indexes(
  thP,
  th,
  replicates = 3,
  pvalue = 0.01,
  states,
  alpha = 0.01,
  p_adjust = "none",
  min_timepoints = 1
)
```

## Arguments

- thP:

  output of output_tc(percent = TRUE). Raw data for procent deuteration
  for time courses

- th:

  output of output_tc(). Raw data for uptake deuteration for time
  courses

- replicates:

  number of replicates in sample. Default set to 3.

- pvalue:

  p-value cutoff. Default set up to 0.01

- states:

  Protein states from the set. As default all states are chosen.

- alpha:

  critical interval, to have more restrictive use lower values,
  default=0.01

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- min_timepoints:

  a peptide is selected when it is significant at at least this many
  time points. Default `1` (any time point, as before).

## Value

Returns indexes of significant peptides: row numbers in the states\[1\]
vs states\[2\] table used by the robot plots. A significant peptide that
is not in that table (e.g. its sequence differs in a mutant) is left
out;
[`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
lists every significant peptide.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
tm_df<-output_tc(filepath=file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written

 # more restictive peptide selection
robot_indexes(thP = tmP_df, th=tm_df, pvalue=0.001, alpha=0.01)
#>  [1]  1  2  3  4  7  8  9 10 11 12 13 15 16 17 18 19 21 22 23 24 26 28
```
