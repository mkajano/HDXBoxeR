# Raw all_data input

The output is being used later for the analysis to compare different
protein states. no initial matching of the files is made, raw data
imported

## Usage

``` r
all_data_long_input(filepath)
```

## Arguments

- filepath:

  filepath to the input file. Input file is All_results table from
  HDX_Examiner, where all the fields are marked for export. Can also be
  a DynamX cluster/state export, or a table from
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  /
  [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md).

## Value

data frame with reorganized in a long format
