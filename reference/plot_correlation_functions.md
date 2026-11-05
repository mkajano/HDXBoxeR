# Plot correlated states

Correlation of two states, the st.dev are shown as ellipsoids.

## Usage

``` r
plot_correlation_functions(df, replicates = 3, ranges, mar_x = 1)
```

## Arguments

- df:

  standard deviation from one sample

- replicates:

  number of replicates. Default set to 3.

- ranges:

  ylim range and and how often ticks

- mar_x:

  x margin variable to be changed if needed

## Value

plot for the correlation

## Examples

``` r
sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
CI_single(s1=sd1, replicates=3)
#> [1] 0.2901757
```
