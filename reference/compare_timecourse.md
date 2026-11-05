# Compare two states over a whole time course

Combines the labeling times of each peptide into whole-time-course
statistics, so a peptide is judged on all its time points together
rather than one time point at a time.

## Usage

``` r
compare_timecourse(
  filepath,
  state1,
  state2,
  times,
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  min_timepoints = 2,
  percent = FALSE,
  replicates = NULL
)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- state1, state2:

  the two protein states (differences are `state2 - state1`).

- times:

  time points to combine. Default: all common time points.

- alpha:

  alpha of the critical intervals. Default `0.01`.

- pvalue:

  p-value cutoff. Default `0.01`.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default) or e.g. `"BH"`; across peptides for the summed test, and
  across peptides within each time for the per-time calls.

- min_timepoints:

  number of time points at which a peptide must pass the per-time rule
  for `sig_min_timepoints`. Default `2`.

- percent:

  value scale as in
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md):
  `FALSE` (Da, default), `TRUE`, `"FD"` or `"per_residue"`.

- replicates:

  optional maximum number of replicates per state.

## Value

data frame with one row per peptide: `Start`, `End`, `Sequence`,
`Charge`, `n_times`, `summed_delta`, `se`, `df`, `p_value`, `p_adj`,
`ci_sum`, `sig_summed`, `n_timepoints_sig`, `sig_min_timepoints`, and
one `sig_<time>` column per time point.

## Details

For each peptide (present at every requested time in both states):

- Summed difference:

  `summed_delta` = sum over time points of (mean state2 - mean state1),
  in Da (or percent). Its standard error combines the per-time standard
  errors, `p_value` is a Welch-type t-test of the sum
  (Welch-Satterthwaite degrees of freedom), and `ci_sum` is a global
  critical interval for summed differences: the per-time pooled Hageman
  & Weis standard errors combined over time points, times
  `|qt(alpha/2, df)|`. `sig_summed` uses the package-wide rule:
  (adjusted) p \< `pvalue` AND \|summed_delta\| \> `ci_sum`.

- Per time point:

  the usual peptide rule at each time (Welch p, adjusted across peptides
  within that time with `p_adjust`, below `pvalue` AND \|delta\| \> that
  time's pooled CI); `n_timepoints_sig` is how many time points pass,
  and `sig_min_timepoints` whether that is at least `min_timepoints`.

## See also

[`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md),
[`woods_CI_plot()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot.md),
[`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
tc <- compare_timecourse(file_nm, "Unbound", "bound")
head(tc[order(tc$p_value), ])
#>    Start End              Sequence Charge n_times summed_delta         se
#> 12    98 115    YSHMVRPGDLTGHSDFHL      2       4    -5.240667 0.10867331
#> 11    95 115 WRFYSHMVRPGDLTGHSDFHL      3       4    -4.834667 0.11218932
#> 17   136 147          IIRLRKGLASRC      3       4    -3.690667 0.10995706
#> 9     89  97             ASVEQFWRF      2       4    -3.061000 0.05518857
#> 21   169 177             AVVSVRFQE      2       4    -1.020333 0.04050789
#> 7     78  88           YEQNIKQIGTF      2       4    -1.439000 0.04927248
#>           df      p_value        p_adj    ci_sum sig_summed n_timepoints_sig
#> 12 11.008622 3.667852e-14 3.667852e-14 0.3434353       TRUE                4
#> 11  8.877092 1.280392e-11 1.280392e-11 0.3434353       TRUE                4
#> 17  9.962907 1.398322e-11 1.398322e-11 0.3434353       TRUE                3
#> 9   7.922340 1.509051e-11 1.509051e-11 0.3434353       TRUE                3
#> 21 11.539764 1.900268e-11 1.900268e-11 0.3434353       TRUE                2
#> 7  10.552846 1.945552e-11 1.945552e-11 0.3434353       TRUE                1
#>    sig_min_timepoints sig_3.00s sig_60.00s sig_1800.00s sig_72000.00s
#> 12               TRUE      TRUE       TRUE         TRUE          TRUE
#> 11               TRUE      TRUE       TRUE         TRUE          TRUE
#> 17               TRUE      TRUE       TRUE        FALSE          TRUE
#> 9                TRUE      TRUE       TRUE         TRUE         FALSE
#> 21               TRUE     FALSE      FALSE         TRUE          TRUE
#> 7               FALSE     FALSE       TRUE        FALSE         FALSE
```
