# Read an HDXExaminer "All results" table

Reads the csv exported by HDXExaminer (Tools \> Export \> All results
table), including exports whose data rows carry one more field than the
header (the "Confidence column shift").

## Usage

``` r
read_hdxexaminer(path)
```

## Arguments

- path:

  path to the csv file.

## Value

data frame of class `hdx_data`.

## See also

[`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
d <- read_hdxexaminer(file_nm)
```
