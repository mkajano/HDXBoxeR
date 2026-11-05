# Global confidence interval treshold from experimental standard deviation for 2 samples.

Calculation of global confidence interval using approach by: Reliable
Identification of Significant Differences in Differential Hydrogen
Exchange-Mass Spectrometry Measurements Using a Hybrid Significance
Testing Approach Tyler S. Hageman and David D. Weis Analytical Chemistry
2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325 The default
alpha is 0.01.

## Usage

``` r
CI_2pts(s1, s2, replicates = 3, alpha = 0.01, n1 = NULL, n2 = NULL)
```

## Arguments

- s1:

  standard deviation from one sample

- s2:

  standard deviation from seconda sample

- replicates:

  number of replicates. A single number, or two numbers (state 1,
  state 2) when the states have different replicate counts. Default set
  to 3.

- alpha:

  probability to calculate t-value in Welch test. Default 0.01

- n1, n2:

  optional per-peptide numbers of replicates behind `s1` and `s2` (e.g.
  from attribute `"n"` of
  [`sd_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timepoint.md)).
  Default: every peptide has `replicates` values.

## Value

treshold for determining significance.

## Details

The t-value is two-tailed, `|qt(alpha/2, df)|`, as in the published
version of the package.

Unequal or missing replicates: each state's pooled SD weights every
peptide's SD by its degrees of freedom (n_i - 1), the standard error of
the difference is `sqrt(sp1^2/n1 + sp2^2/n2)` with n1, n2 the average
number of replicates per peptide in each state, and the t-value uses df
= n1 + n2 - 2. With equal, complete replicates this is exactly the
Hageman & Weis formula.

## Examples

``` r
sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
sd2<-data.frame(c(0.18, 0.11, 0.13, 0.08, 0.11, 0.06))
CI_2pts(s1=sd1, s2=sd2, replicates=3, alpha=0.01)
#> [1] 0.4272419
```
