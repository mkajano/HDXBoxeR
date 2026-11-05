# Calculate what would be the standard deviation based on the CI

Calculate the CI to standard devaation, that can be used for the
ellipsoid function Calculate standard deviation from the CI using the
equation as shown in:
SD=sqrt(replicates)*(uper_limit-lower_limit)/(2*t), t=qt(alpha/2, df)
for the two-sided (1-alpha) CI returned by t.test()

## Usage

``` r
sd_ci_timepoint(df, replicates = 3, alpha = 0.01)
```

## Arguments

- df:

  input from tp function

- replicates:

  number of replicates in sample. Default set to 3.

- alpha:

  what alpha (CI cutoff should be used)

## Value

data frame with standard deviation calculated from the CI
