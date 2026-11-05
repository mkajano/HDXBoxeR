# Returns dataframe with peptides which exhibit significant difference between two sets

Function to help decide which peptides will be drawn on Robot plots.

## Usage

``` r
robot_indexes_df(
  thP,
  th,
  replicates = 3,
  pvalue = 0.01,
  alpha = 0.01,
  states,
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

- alpha:

  cutoff for Critical interval. Default=0.01

- states:

  Protein states from the set. As default all states are chosen.

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- min_timepoints:

  a peptide is selected when it is significant at at least this many
  time points. Default `1` (any time point, as before); `2` or more
  ignores one-off calls at a single time point.

## Value

Returns dataframe listing peptides that are significantly different
between sets.

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
robot_indexes_df(thP = tmP_df, th=tm_df, pvalue=0.001, alpha=0.01)
#>    nb Start End              Sequence Charge
#> 1   1    40  58   GPLGSKAVVPGPAEHPLQY      2
#> 2   2    40  59  GPLGSKAVVPGPAEHPLQYN      2
#> 3   3    47  63     VVPGPAEHPLQYNYTFW      2
#> 4   4    62  77      FWYSRRTPGRPTSSQS      3
#> 5   5    78  88           YEQNIKQIGTF      2
#> 6   6    82  88               IKQIGTF      2
#> 7   7    89  97             ASVEQFWRF      2
#> 8   8    95 107         WRFYSHMVRPGDL      3
#> 9   9    95 115 WRFYSHMVRPGDLTGHSDFHL      3
#> 10 10    98 115    YSHMVRPGDLTGHSDFHL      2
#> 11 11   114 124           HLFKEGIKPMW      2
#> 12 12   116 124             FKEGIKPMW      2
#> 13 13   116 135  FKEGIKPMWEDDANKNGGKW      3
#> 14 14   136 147          IIRLRKGLASRC      3
#> 15 15   136 149        IIRLRKGLASRCWE      3
#> 16 16   140 147              RKGLASRC      2
#> 17 17   169 177             AVVSVRFQE      2
#> 18 18   176 183              QEDIISIW      2
#> 19 19   176 193    QEDIISIWNKTASDQATT      2
#> 20 20   178 193      DIISIWNKTASDQATT      2
#> 21 21   183 193           WNKTASDQATT      2
#> 22 22   194 200               ARIRDTL      2
```
