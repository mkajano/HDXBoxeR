# Returns data frame with difference of averages between State1 and other states provided.

Returns average difference data.frame. Every state is compared to the
first state in `df` (state - first state). To use another reference, put
it first with `output_tp(states = ...)`.

## Usage

``` r
dif_ave(df)
```

## Arguments

- df:

  output of ave_timepoint() on an output_tp() table.

## Value

Data.frame with difference values btw control and other protein states.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
d<-dif_ave(ave_timepoint(a))
## bound as the reference instead
d2<-dif_ave(ave_timepoint(output_tp(file_nm, states=c("bound", "Unbound"))))
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: bound=3 Unbound=3
#> no csv file written
```
