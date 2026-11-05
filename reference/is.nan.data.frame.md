# Checks for NaN is data.frame

Function by Hong Ooi;
https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame

## Usage

``` r
# S3 method for class 'data.frame'
is.nan(x)
```

## Arguments

- x:

  Data frame to be checked for NaN

## Value

logical. Returns info if data.frame contains NaNs.

## Examples

``` r
## this function will overwrite the is.nan function that works only on vectors and matrices
df<-data.frame(c(0,NaN), c(1, 2))
is.nan(df)
#>      c.0..NaN. c.1..2.
#> [1,]     FALSE   FALSE
#> [2,]      TRUE   FALSE
df[is.nan(df)]<- 0
```
