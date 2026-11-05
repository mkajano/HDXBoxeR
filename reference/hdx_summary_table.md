# Summary table for HDX-MS reporting (Masson et al. 2019)

Builds the summary table recommended by Masson et al. (2019, Nat Methods
16:595) for publishing HDX-MS data: one row per item, one column per
protein state. Items computed from the data: time course, controls,
replicates, number of peptides, sequence coverage, average peptide
length, redundancy, back-exchange (mean and range), repeatability
(pooled replicate SD, Da) and the significance threshold: the critical
interval used when each state is compared with the first state (the same
values as in
[`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md)).
Items that only you know (labeling buffer, pH, temperature, quench,
protease, LC/MS, software ...) are added from `metadata`.

## Usage

``` r
hdx_summary_table(
  filepath,
  metadata = list(),
  replicates = NULL,
  Dfact = 0.85,
  alpha = 0.01,
  pvalue = 0.01,
  protein_range = NULL,
  file = NULL,
  states = NULL
)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- metadata:

  named list of extra items, e.g.
  `list("HDX reaction details" = "20 mM HEPES pD 7.4, 90% D2O, 20 C")`.
  A single value applies to every state; a vector gives one value per
  state. They are listed first.

- replicates:

  optional maximum number of replicates per state. Default `NULL` uses
  every replicate.

- Dfact:

  D2O fraction of the labeling buffer for back-exchange. Default `0.85`.

- alpha:

  alpha of the critical interval. Default `0.01`.

- pvalue:

  p-value cutoff, reported with the threshold. Default `0.01`.

- protein_range:

  residue range for the coverage, as in
  [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md).
  Default: first to last covered residue.

- file:

  optional output file: `.csv` or `.md` (markdown table).

- states:

  optional states to include (the first one is the reference). Default:
  every state.

## Value

data frame with columns `Item` and one per state (invisibly when `file`
is given).

## See also

[`all_summary()`](https://mkajano.github.io/HDXBoxeR/reference/all_summary.md),
[`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md),
[`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
hdx_summary_table(file_nm, metadata = list(
  "HDX reaction details" = "PBS pD 7.4, 85% D2O, 22 C",
  "Data analysis software" = "HDExaminer 3; HDXBoxeR"))
#>                               Item                                  Unbound
#> 1             HDX reaction details                PBS pD 7.4, 85% D2O, 22 C
#> 2           Data analysis software                   HDExaminer 3; HDXBoxeR
#> 3                  HDX time course          3.00s 60.00s 1800.00s 72000.00s
#> 4                     HDX controls undeuterated (0s); fully deuterated (FD)
#> 5             Number of replicates                                        3
#> 6               Number of peptides                                       37
#> 7            Sequence coverage (%)                  94.09 (residues 40-225)
#> 8           Average peptide length                                    13.86
#> 9                       Redundancy                                     2.93
#> 10 Back-exchange (%): mean (range)                     27.01 (1.31 - 51.48)
#> 11   Repeatability (pooled SD, Da)                                    0.067
#> 12  Significant differences in HDX                          reference state
#>                                                                               bound
#> 1                                                         PBS pD 7.4, 85% D2O, 22 C
#> 2                                                            HDExaminer 3; HDXBoxeR
#> 3                                                   3.00s 60.00s 1800.00s 72000.00s
#> 4                                          undeuterated (0s); fully deuterated (FD)
#> 5                                                                                 3
#> 6                                                                                34
#> 7                                                           94.09 (residues 40-225)
#> 8                                                                             13.29
#> 9                                                                              2.58
#> 10                                                            35.44 (11.73 - 99.88)
#> 11                                                                            0.076
#> 12 |difference| > 0.20-0.35 Da (critical interval, alpha = 0.01) and Welch p < 0.01
```
