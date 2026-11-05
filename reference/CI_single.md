# Global confidence interval treshold from experimental standard deviation for 1 sample

Calculation of global confidence interval using approach by: Reliable
Identification of Significant Differences in Differential Hydrogen
Exchange-Mass Spectrometry Measurements Using a Hybrid Significance
Testing Approach Tyler S. Hageman and David D. Weis Analytical Chemistry
2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325 The default
alpha is 0.01.

## Usage

``` r
CI_single(s1, replicates = 3, alpha = 0.01)
```

## Arguments

- s1:

  standard deviation from one sample

- replicates:

  number of replicates. Default set to 3.

- alpha:

  probability to calculate t-value in Welch test. Default 0.01

## Value

treshold for determining significance.

## Examples

``` r
sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
CI_single(s1=sd1, replicates=3, alpha=0.01)
#> [1] 0.2901757
```
