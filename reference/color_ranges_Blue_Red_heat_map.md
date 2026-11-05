# Returns color pallete from red to blue with number of colors for defined ranges

Returns color pallete from red to blue with number of colors for defined
ranges

## Usage

``` r
color_ranges_Blue_Red_heat_map(ranges, col, colors_initial)
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
color_ranges_Blue_Red_heat_map(ranges=c(-Inf, -100, -50, 0, 50, 100, Inf), col="white")
#> [1] "white"   "#A50F15" "#EF3B2C" "#FC9272" "#9ECAE1" "#4292C6" "#08519C"
```
