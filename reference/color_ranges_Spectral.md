# Returns Spectral pallette with colors matching defined ranges

Spectral pallette for timecourse data

## Usage

``` r
color_ranges_Spectral(ranges, col, colors_initial)
```

## Arguments

- ranges:

  vector of numbers. Should have the same mumber of positive and
  negative values and contain 0.

- col:

  additional color that should be first in the pallette.

- colors_initial:

  deprecated, use `col`.

## Value

color scheme for number

## Examples

``` r
color_ranges_Spectral(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), col="white")
#> [1] "white"   "#5E4FA2" "#66C2A5" "#E6F598" "#FEE08B" "#F46D43" "#9E0142"
```
