# Hidden function from qpcR package, typical usage as qpcR:::cbind.na

Combine data of unequal row length avoiding repetition or errors by
filling with NAs. In contrast to classical cbind, cbind.na can be used
to combine data such as

## Usage

``` r
qpcr.cbind.na(..., deparse.level = 1)
```

## Arguments

- ...:

  vectors

- deparse.level:

  set to 1 as default

## Value

data frame with NA

## Examples

``` r
qpcr.cbind.na(1:10, 1:3)
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    2    2
#>  [3,]    3    3
#>  [4,]    4   NA
#>  [5,]    5   NA
#>  [6,]    6   NA
#>  [7,]    7   NA
#>  [8,]    8   NA
#>  [9,]    9   NA
#> [10,]   10   NA
```
