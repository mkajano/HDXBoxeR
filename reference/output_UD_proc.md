# Prepares output for HDX-MS Undeuterated data for procent deuteration.

Returns the in-exchange control (`Deut.Time` 0, e.g. `"0s"`) as percent
deuteration: one row per peptide and protein state, one column per
replicate.

## Usage

``` r
output_UD_proc(filepath)
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

data frame with reorganized data where in columns is procent deuteration
for Protein States.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_UD_proc(file_nm)
```
