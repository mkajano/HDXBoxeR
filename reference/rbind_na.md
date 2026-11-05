# bind non equal row

kmezhoud/canceR: A Graphical User Interface for accessing and modeling
the Cancer Genomics Data of MSKCC
https://rdrr.io/github/kmezhoud/canceR/src/R/rbind.na.R

## Usage

``` r
rbind_na(..., deparse.level = 1)
```

## Arguments

- ...:

  (generalized) vectors or matrices.

- deparse.level:

  integer controlling the construction of labels in the case of
  non-matrix-like arguments (for the default method): deparse.level = 0
  constructs no labels; the default, deparse.level = 1 or 2 constructs
  labels from the argument names.

## Value

a data frame with merged rows

## Examples

``` r
row1 <- c("a","b","c","d")
row2 <- c("A", "B", "C")
row3 <- rbind_na(row1, row2)
```
