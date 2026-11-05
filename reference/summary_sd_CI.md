# Provides summary table with Critical interval and standard deviation within the set.

For each protein state: the pooled replicate standard deviation of all
peptides and time points (each peptide's SD weighted by its replicates
minus one), and the critical interval of that state alone
([`CI_single()`](https://mkajano.github.io/HDXBoxeR/reference/CI_single.md))
using the state's average number of replicates.

## Usage

``` r
summary_sd_CI(filepath, replicates = NULL, alpha = 0.01)
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

- alpha:

  critical interval, to have more restrictive use lower values,
  default=0.01

## Value

Returns summary table.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- summary_sd_CI(file_nm, alpha=0.01)
```
