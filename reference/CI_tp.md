# Global confidence interval treshold from experimental standard deviation

Calculation of global confidence interval using approach by for all
protein states compared to first state in the data.frame. Reliable
Identification of Significant Differences in Differential Hydrogen
Exchange-Mass Spectrometry Measurements Using a Hybrid Significance
Testing Approach Tyler S. Hageman and David D. Weis Analytical Chemistry
2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325

## Usage

``` r
CI_tp(df, replicates = 3, alpha = 0.01)
```

## Arguments

- df:

  standard deviation dataframe.

- replicates:

  number of replicates. Default set to 3.

- alpha:

  significance level. Set as default to 0.01

## Value

treshold for determining significance.

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
CI_tp(df=sd, replicates=3, alpha=0.01 )
#> [1] 0.2706825
CI_tp(sd)
#> [1] 0.2706825
```
