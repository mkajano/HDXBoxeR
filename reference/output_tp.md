# Prepares output for HDX-MS for the deuteration uptake or percent deuteration for the time points.

Returns a data frame organized for additional analysis. In columns are
deuteration uptake or percent deuteration data for the given protein
states. Function allows for writing csv with data, matching sequences of
peptide. Protein.States, Deut.times, or number of replicates can be
specified.

## Usage

``` r
output_tp(
  filepath,
  replicates,
  states,
  times,
  seq_match = FALSE,
  csv = "NA",
  percent = FALSE,
  min_replicates = 2,
  .arg_df = NULL
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

  optional maximum number of replicates per state and time. By default
  every replicate is used, so states may have different numbers of
  replicates (e.g. 3 and 4); missing measurements are `NA`. Replicates
  are numbered in the order their experiments appear in the file.

- states:

  function allows to choose what states should be used for analysis.
  Default all states are used.

- times:

  lists the deuteration times to be used in analysis. Default all states
  used.

- seq_match:

  How peptides are matched between states. `FALSE` (default): a peptide
  is kept only if the same Sequence, Start, End and Charge are found in
  every state (exact sequence identity). `TRUE`: peptides are matched by
  position only (Start, End and Charge), ignoring Sequence, so peptides
  covering e.g. a point mutation are kept; where the sequences differ
  between states they are joined with "." in the `Sequence` column.

- csv:

  Flag allowing saving the output as csv. With default csv="NA", data is
  not saved. If csv output is desided, provide output name.

- percent:

  What the values are. `FALSE` (default): deuterium uptake in Da.
  `TRUE`: percent deuteration as reported by the input (HDXExaminer's
  `Deut %`, i.e. relative to the theoretical maximum; for DynamX 100 x
  uptake / `MaxUptake`). `"FD"`: percent of the measured fully
  deuterated control, 100 x uptake / FD uptake of the same peptide and
  state (FD replicates averaged, charge states pooled); this corrects
  for back-exchange and D2O content. Peptides without a usable FD value
  (missing, or FD uptake \<= 0) keep their row with empty (`NA`) values.
  `"per_residue"`: uptake per exchangeable amide, Da /
  [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
  (length - 2 - prolines after the first residue), so peptide length
  does not affect the critical-interval half of the significance rule
  (`TRUE` and `"FD"` are length-normalized too). The p-values are the
  same on every scale; only the critical interval and the differences
  change.

- min_replicates:

  a peptide is kept only if it has at least this many replicate values
  in every state at that time. Default `2`, the minimum for a standard
  deviation and a t-test.

- .arg_df:

  advanced/internal: a pre-computed
  [`arg_df()`](https://mkajano.github.io/HDXBoxeR/reference/arg_df.md)
  data.frame for this same `filepath`, to skip re-reading and re-parsing
  the file from disk. Default `NULL` reads `filepath` as usual - only
  pass this when calling `output_tp()` many times for the same file
  (e.g. once per pair of states) and `arg_df(filepath)` has already been
  computed once by the caller; passing a data.frame that doesn't match
  `filepath` will silently produce wrong results, since nothing here
  re-checks it.

## Value

data frame with reorganized data where in columns is the deuteration
uptake for Protein States. Attributes `"replicates"` (number of
replicate columns per state), `"states"`, `"times"` and `"percent"`
record how it was built.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(filepath=file_nm) ###all default parameters used
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written


# all possible flags listed & percent deuteration output,
# with sequences matching for protein states.

a<-output_tp(filepath=file_nm, replicates=3, states=c("bound", "Unbound"),
times=c("3.00s", "72000.00s"), seq_match=TRUE, csv="NA", percent=TRUE)
#> 14 peptide row(s) dropped: 14 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: bound=3 Unbound=3
#> no csv file written
```
