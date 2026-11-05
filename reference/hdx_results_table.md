# Tidy table of peptide-level results

One row per peptide, time point and compared pair of states, with the
means, SDs, replicate counts, difference, Welch p-value, adjusted
p-value, critical interval and significance call - the package-wide rule
(p \< `pvalue` AND \|difference\| \> pooled CI). Suitable for
supplementary tables and for further analysis outside HDXBoxeR. The
analysis parameters, package version and date are attached as attribute
`"hdxboxer"`; with `file`, the table is written as csv together with a
`<file>_parameters.txt` holding the parameters and
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

## Usage

``` r
hdx_results_table(
  filepath,
  states,
  times,
  reference,
  all_pairs = FALSE,
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  percent = FALSE,
  replicates = NULL,
  file = NULL
)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- states:

  states to compare. Default: all.

- times:

  time points. Default: all.

- reference:

  reference state; every other state is compared with it (difference =
  other - reference). Default: the first state. `reference = NULL` or
  `all_pairs = TRUE` compares every pair.

- all_pairs:

  `TRUE` compares every pair of states instead.

- alpha:

  alpha of the pooled critical interval. Default `0.01`.

- pvalue:

  p-value cutoff. Default `0.01`.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default) or e.g. `"BH"`, across the peptides of each pair of states
  at each time point.

- percent:

  value scale as in
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md):
  `FALSE` (Da, default), `TRUE`, `"FD"` or `"per_residue"`.

- replicates:

  optional maximum number of replicates per state.

- file:

  optional csv file to write.

## Value

data frame (invisibly when `file` is given).

## See also

[`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md),
[`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
res <- hdx_results_table(file_nm)
head(res[res$is_sig, ])
#>    Deut.Time Start End             Sequence Charge  state1 state2    mean1
#> 1      3.00s    40  58  GPLGSKAVVPGPAEHPLQY      2 Unbound  bound 8.327667
#> 2      3.00s    40  59 GPLGSKAVVPGPAEHPLQYN      2 Unbound  bound 8.717333
#> 3      3.00s    47  63    VVPGPAEHPLQYNYTFW      2 Unbound  bound 3.974667
#> 8      3.00s    82  88              IKQIGTF      2 Unbound  bound 1.756000
#> 9      3.00s    89  97            ASVEQFWRF      2 Unbound  bound 1.396000
#> 10     3.00s    95 107        WRFYSHMVRPGDL      3 Unbound  bound 3.610333
#>           sd1 n1    mean2        sd2 n2      delta      p_value        p_adj
#> 1  0.13838473  3 7.013667 0.05379901  3 -1.3140000 1.298206e-03 1.298206e-03
#> 2  0.15430597  3 7.468667 0.10442382  3 -1.2486667 6.305205e-04 6.305205e-04
#> 3  0.10868456  3 2.785667 0.04764801  3 -1.1890000 6.958101e-04 6.958101e-04
#> 8  0.07158212  3 1.264667 0.05953430  3 -0.4913333 9.252337e-04 9.252337e-04
#> 9  0.04167733  3 0.827000 0.06199193  3 -0.5690000 4.138932e-04 4.138932e-04
#> 10 0.05381759  3 1.558667 0.08418036  3 -2.0516667 1.706338e-05 1.706338e-05
#>          ci is_sig
#> 1  0.350034   TRUE
#> 2  0.350034   TRUE
#> 3  0.350034   TRUE
#> 8  0.350034   TRUE
#> 9  0.350034   TRUE
#> 10 0.350034   TRUE
attr(res, "hdxboxer")$parameters$p_adjust
#> [1] "none"
```
