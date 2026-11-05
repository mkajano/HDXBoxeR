# Prepares the long input of the input data

The output is being used later for the analysis to compare different
protein states. input parameters match the output_tp parameters

## Usage

``` r
boxer_input_to_long(
  filepath,
  replicates,
  states,
  times,
  seq_match = FALSE,
  csv = "NA",
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

- replicates:

  number of replicates to be used in analysis. The function takes number
  of replicates up to specified number. If no argument provided number
  maximal common number of replicates it used.

- states:

  function allows to choose what states should be used for analysis.
  Default all states are used.

- times:

  lists the deuteration times to be used in analysis. Default all states
  used.

- seq_match:

  Flag allows to choose if the peptide sequences should be matched
  between states. seq_match=FALSE signifies no sequence matching,
  seq_match=TRUE states that the sequences are matched between the sets.

- csv:

  Flag allowing saving the output as csv. With default csv="NA", data is
  not saved. If csv output is desided, provide output name.

- percent:

  Flag allowing to choose output as deteuration uptake (FALSE) or
  percent deuteration (TRUE). Default deuteration uptake.

## Value

data frame with reorganized in a long format
