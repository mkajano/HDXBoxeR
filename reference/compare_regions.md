# Region-level ANOVA: compare protein states within residue regions

Tests whether protein states differ in deuterium uptake within each
region you define. Each peptide acts as its own control: the state
effect is a two-way ANOVA `uptake ~ peptide + state` on the per-peptide
mean uptake (replicates averaged), so it is the uptake difference
*within* peptides, averaged over the region's peptides, judged against
how consistent that difference is from peptide to peptide (peptides are
the replication unit, as in a randomized-block / repeated-measures
design; with two states this is a paired t-test of the per-peptide
differences). This is far more powerful than pooling peptide means into
one distribution per state (as
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
must do when peptide maps differ), because differences between peptides
(length, location) are removed before the states are compared.

## Usage

``` r
compare_regions(
  filepath,
  times,
  regions,
  states = NULL,
  reference = NULL,
  by = c("sequence", "position"),
  replicates = NULL,
  percent = FALSE,
  time_model = FALSE,
  p_adjust = "none"
)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- times:

  one or more `Deut.Time` values. Default: all.

- regions:

  named list of `c(start, end)` residue ranges, e.g.
  `list(NTR = c(1, 85), ACD = c(86, 165))`. A peptide belongs to a
  region when it lies entirely inside it.

- states:

  states to compare (two or more). Default: all.

- reference:

  optional state; pairwise contrasts are then only `reference` vs each
  other state.

- by:

  peptide matching between states: `"sequence"` (default) or
  `"position"`, as in
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md).

- replicates:

  optional maximum number of replicates per state.

- percent:

  value scale, as in
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md):
  `FALSE` (Da, default), `TRUE`, `"FD"` or `"per_residue"`.

- time_model:

  `FALSE` (default): one model per region and time point. `TRUE`: one
  model per region over all `times` (needs two or more).

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default) or e.g. `"BH"`.

## Value

a list of two data frames. `summary`: one row per region (and time
point): `region`, `Deut.Time` (`"all"` with `time_model = TRUE`),
`n_peptides`, `n_obs`, `F_state`, `df_state`, `df_resid`, `p_state`,
`p_state_adj`, `p_state_x_peptide`, and with `time_model = TRUE` also
`p_state_x_time`. `pairs`: one row per region, time and state pair:
`estimate` (Da, or percent, state2 - state1), `se`, `p_value`, `p_adj`.

## Details

Requires the same peptides in every compared state (matched peptides, as
for
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md));
peptides missing from any state are left out of that region. For states
with different peptide maps use
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md).

For each region (and time point) the result reports:

- an F-test of the state effect (do the states differ at all in this
  region; with two states it is the test of their difference),

- an F-test of the state x peptide interaction on the replicate values
  (does the state effect differ between peptides of the region by more
  than replicate noise - a significant value means the effect is
  localized or mixed in direction, so the average shift alone does not
  describe the region; look at the peptides with
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
  or
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)),

- pairwise state contrasts: average shift (Da, state2 - state1) per
  peptide with its standard error and t-test p-value from the model.

With several time points and `time_model = TRUE`, one model per region
is fitted over all of them on the per-peptide-and-time means,
`uptake ~ peptide * time + state`, testing the state effect averaged
over time, plus whether it changes with time (state x time).

A region needs at least two peptides. Every peptide weighs equally,
whatever its replicate noise. p-values are raw by default; `p_adjust`
corrects across all regions (and time points), and separately across all
pairwise contrasts.

## See also

[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
(peptide by peptide),
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
(different peptide maps),
[`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)
(regions found from the data)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
res <- compare_regions(file_nm, times = "1800.00s",
                       regions = list(nterm = c(1, 100), cterm = c(101, 225)))
res$summary
#>   region Deut.Time n_peptides n_obs   F_state df_state df_resid    p_state
#> 1  nterm  1800.00s          9    54  4.012432        1        8 0.08012977
#> 2  cterm  1800.00s         20   120 17.078329        1       19 0.00056602
#>   p_state_adj p_state_x_peptide
#> 1  0.08012977      7.621526e-28
#> 2  0.00056602      2.886840e-24
res$pairs
#>   region Deut.Time  state1 state2   estimate         se    p_value      p_adj
#> 1  nterm  1800.00s Unbound  bound -0.3596296 0.17953604 0.08012977 0.08012977
#> 2  cterm  1800.00s Unbound  bound -0.2543000 0.06153521 0.00056602 0.00056602
```
