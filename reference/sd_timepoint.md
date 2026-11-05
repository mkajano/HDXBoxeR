# Returns standard deviation for dataframe.

Calculates the standard deviation of the replicates of every state (or
time, for
[`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md)
tables) for each peptide. Replicate blocks are found from the column
names, so states may have different numbers of replicates and missing
replicates (`NA`) are ignored; a peptide with fewer than two values in a
block gets `NA`.

## Usage

``` r
sd_timepoint(df, replicates = 3)
```

## Arguments

- df:

  output from functions output_tp or output_tc.

- replicates:

  only used for tables whose value columns are not named
  `<state>_..._<replicate>`: the number of replicates per block. Default
  3.

## Value

Data.frame with standard deviation. The per-peptide number of replicates
behind each SD is attached as attribute `"n"`.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
sd<-sd_timepoint(df=a, replicates=3)
```
