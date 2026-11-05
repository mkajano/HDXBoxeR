# Number of exchangeable protons

Provides a vector with number of exchangeable protons, calculated from
the input table. Number of protons calculated as peptide_length - 2 -
number of Prolines in the peptide that are not in the first position.
Note: if the peptide's *second* residue is a Proline, it is excluded
both by the generic "-2" and by the Proline count, so the result is one
lower than `peptide_exch_residues()`'s residue count for that peptide.

## Usage

``` r
nb_exch_deut(df)
```

## Arguments

- df:

  standard deviation from one sample

## Value

vector with number of exchangeable protons

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
nb_exch_deut(a)
#>   [1] 13 14 12 12 11 10  9  5  7 10 18 15  8 19  6 17 10 12  6  6  7  6 16 14 10
#>  [26]  9 16  5  9  6 10  5 13 14 12 12 11 10  9  5  7 10 18 15  8 19  6 17 10 12
#>  [51]  6  6  7  6 16 14 10  9 16  5  9  6 10  5 13 14 12 12 11 10  9  5  7 10 18
#>  [76] 15  8 19  6 17 10 12  6  6  7  6 16 14 10  9 16  5  9  6 10  5 13 14 12 12
#> [101] 11 10  9  5  7 10 18 15  8 19  6 17 10 12  6  6  7  6 16 14 10  9 16  5  9
#> [126]  6 10  5
```
