# Prepares output for HDX-MS Undeuterated sample data.

Returns the in-exchange control (`Deut.Time` 0, e.g. `"0s"` or
`"0.00s"`): one row per peptide and protein state, one column per
replicate (uptake in Da). Missing uptake (HDXExaminer writes `n/a` for
the undeuterated sample) is set to 0, since the undeuterated sample is
the reference.

## Usage

``` r
output_UD(filepath)
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

data frame with reorganized data where in columns is uptake data for
Protein States.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_UD(file_nm)
```
