# Allows for selecting some peptide from input data

Function allows for picking indices from the inputs based on: peptide
start or end residue, length, state or timepoint. If parameters set to
NA, condition is skipped.

## Usage

``` r
select_indices(df, start = NA, end = NA, length = NA, times = NA, states = NA)
```

## Arguments

- df:

  input file (output of output_tc or output_tp)

- start:

  provide number for the staring residue, default NA

- end:

  provide number for the end residue, default NA

- length:

  provide max length of the peptide

- times:

  timepoints, only for the output_tp functions

- states:

  states, only for the output_tc functions

## Value

Row indices of the peptides that are fulfilling the conditions required.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
indb<-select_indices(a,length=12, start=100, end=200)
smaller_df<-a[indb,]
```
