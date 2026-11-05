# Preparatory function writing pymol scripts

Function rearrange vector to string by adding + sign between the
numbers.

## Usage

``` r
pymol_str(ind1)
```

## Arguments

- ind1:

  vector of numbers (residues)

## Value

string with + as a separator.

## Examples

``` r
res<-c(1,5, 19, 100, 109)
pymol_str(res)
#> [1] "1+5+19+100+109"
```
