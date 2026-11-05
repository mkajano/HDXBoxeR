# Critial interval calculation two sets of timecourses

Preparatory function for calculation of pvalue between sets.

## Usage

``` r
CI_tc(sd_c, sd_v, replicates = 3, alpha = 0.01)
```

## Arguments

- sd_c:

  dataframe of control

- sd_v:

  dataframe for variant

- replicates:

  number of replicates. Default set to 3.

- alpha:

  pvalue cutoff. Default set to 0.01

## Value

Critical interval for 2 sets
