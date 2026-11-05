# Prepares output with HDX-MS data for publications

Format prepared based of example from: Masson, G.R., Burke, J.E., Ahn,
N.G. et al. Recommendations for performing, interpreting and reporting
hydrogen deuterium exchange mass spectrometry (HDX-MS) experiments. Nat
Methods 16, 595–602 (2019). https://doi.org/10.1038/s41592-019-0459-y It
generates csv file in format ready for publication of the data.

## Usage

``` r
output_prep(
  filepath,
  output_name,
  states,
  replicates = NULL,
  times,
  percent = FALSE
)
```

## Arguments

- filepath:

  filepath to the input file. Input file is All_results table from
  HDX_Examiner, where all the fields are marked for export. Can also be
  a DynamX cluster/state export, or a table from
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  /
  [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md).

- output_name:

  Name of output file. It has to be csv file

- states:

  function allows to choose what states should be used for analysis.
  Default all states are used.

- replicates:

  optional maximum number of replicates per state and time. Default
  `NULL` uses every replicate; states may differ.

- times:

  lists the deuteration times to be used in analysis. Default all states
  used.

- percent:

  return either uptake or percent deuteration, default=FALSE, return
  uptake

## Value

Writes the table to `output_name` and returns it invisibly.
`Mean.Peptide_Mass_[Da]` is the mean measured neutral peptide mass,
Charge x (centroid m/z - proton mass); earlier versions reported the
centroid m/z under this name.

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
output_prep(filepath=file_nm, output_name=tempfile())
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> Number of replicates used: Unbound=3 bound=3
# }
```
