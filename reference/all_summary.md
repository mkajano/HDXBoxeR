# Returns full summary table.

Returns summary data. Function returns: Protein states, timepoints,
number of replicates, \# peptides, % coveregae, average peptide length
and redundancy. backexchange calculations (average and range), Critical
interval and standard deviation. Function requires undeuterated and
Fully deuterated sets marked in Deut.time as 0s and FD respectively.

## Usage

``` r
all_summary(
  filepath,
  replicates = NULL,
  Dfact = 0.85,
  alpha = 0.01,
  protein_range = NULL
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

  optional maximum number of replicates per state. Default `NULL` uses
  every replicate.

- Dfact:

  Dfact is the fraction of D/H in the labeling buffer used. Default set
  up to 0.85

- alpha:

  Critical interval, level of stringency

- protein_range:

  residue range for the coverage, as in
  [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md).

## Value

Returns summary table. Back-exchange columns are `NA` when the input has
no usable undeuterated/fully deuterated controls.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- all_summary(file_nm, Dfact=0.85)
```
