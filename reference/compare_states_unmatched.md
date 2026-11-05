# Compare pooled peptide-level HDX uptake across states without matching peptides

A deliberately different tool from
[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md):
instead of requiring the same peptide to be identified in every state
being compared, this pools *every* peptide in a chosen residue range
within each state into one distribution, and asks whether that region's
typical exchange level differs between states - no peptide
correspondence required at all. Useful when peptide identification
differs enough between states (chimeras, point mutants,
differently-digested runs) that peptide-matched comparison would discard
too much data, at the cost of a real confound: two states can appear to
differ purely because their peptide *maps* differ (different
lengths/boundaries), not because the underlying exchange does. See
Details.

## Usage

``` r
compare_states_unmatched(
  filepath,
  times,
  regions,
  states = NULL,
  value = c("raw", "percent", "frac_deut"),
  method = c("peptide", "classic", "resample"),
  reference = NULL,
  replicates = 3,
  pvalue = 0.01,
  alpha = 0.01,
  p_adjust = "none",
  n_perm = 1000,
  n_boot = 1000,
  seed = NULL,
  p_adjust_method
)
```

## Arguments

- filepath:

  filepath to the input file (read via
  [`all_data_long_input()`](https://mkajano.github.io/HDXBoxeR/reference/all_data_long_input.md) -
  no peptide matching between states happens at any stage).

- times:

  a single `Deut.Time` value to compare at.

- regions:

  required named list of `c(start, end)` residue ranges, e.g.
  `list(NTR = c(1, 85), ACD = c(86, 165))` - see Details for why
  whole-protein pooling isn't offered.

- states:

  character vector of states to compare. Default `NULL` uses every state
  found in `filepath`.

- value:

  one of `"raw"` (default), `"percent"`, `"frac_deut"` - which column
  drives the comparison. See "What gets compared: `value`".

- method:

  one of `"peptide"` (default), `"classic"`, `"resample"` - see Details.

- reference:

  optional state name - if given, only `reference` vs. each other state
  is tested, instead of every pair.

- replicates:

  not used: every replicate of every state is used. Kept so that older
  calls still work.

- pvalue:

  significance cutoff. This function does not compute `is_sig` for the
  pooled `shift_hl`/`p_value` result itself (callers can filter
  `p_adj < pvalue` there themselves), but `pvalue` does gate the
  per-peptide `frac_peptides_significant` diagnostic - see "Pooled
  significance is not peptide significance". Default `0.01`.

- alpha:

  alpha of the pooled Hageman & Weis critical interval used by the
  per-peptide diagnostic. Default `0.01`.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default) or e.g. `"BH"`, applied across every region-by-pair test
  together, and, for the per-peptide diagnostic, across the peptides of
  each pair of states.

- n_perm:

  number of label permutations for `method = "resample"`. Default
  `1000`.

- n_boot:

  number of bootstrap resamples for `method = "resample"`'s confidence
  interval. Default `1000`.

- seed:

  optional random seed for `method = "resample"`, making the permutation
  p-values and bootstrap interval reproducible. The caller's own
  random-number stream is left unchanged. Default `NULL`.

- p_adjust_method:

  deprecated; earlier versions used `p_adjust_method = "BH"` by default.
  Supplying it sets `p_adjust`, with a warning.

## Value

a list with `value` (the `value` argument actually used, echoed back for
plotting/labeling) and three data.frames: `summary` (one row per
region-by-state-pair: `region`, `state1`, `state2`, `n1`, `n2`,
`shift_hl`, `ci_low`, `ci_high` - `NA` unless `method = "resample"` -
`p_value`, `p_omnibus`, `p_adj`, `n_peptides_tested`,
`n_peptides_significant`, `frac_peptides_significant` - see "Pooled
significance is not peptide significance"), `coverage` (one row per
region-by-state: `n_peptides`, `mean_length`, and
`frac_shared_all_states` - the fraction of that region's peptide
identities, by Start/End/Sequence, present in *every* state being
compared, a diagnostic for how much the matched-vs-unmatched confound
applies here), and `data` (the pooled `deut_value` values the tests
actually ran on - one row per region-by-state-by-peptide, or per
replicate under `method = "classic"` - for plotting with
[`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md)),
and `peptides` (the matched per-peptide results behind the diagnostic:
`delta`, `p_value`, `p_adj`, `ci` and `significant` per region, pair and
peptide).

## This is exploratory, not confirmatory

Unlike
[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md),
results here are `shift_hl`/`ci_low`/ `ci_high`/`p_value` - deliberately
not `is_sig`/`delta` - so this function's output is never mistaken for
the package's confirmatory Hageman & Weis criterion. Always inspect the
returned `coverage` table before trusting a result: if
`frac_shared_all_states` is low for a region, the states' peptide maps
differ enough there that some of any observed shift may reflect peptide
composition, not biology.

## Why `regions` is required

Pooling every peptide in the whole protein mixes a strongly bimodal mix
of protected-core and exposed-loop behavior; the pooled median is then
governed mostly by which peptides happened to be identified (composition
has more sequence available in exposed regions to fully deuterate), not
by biology. Restricting to one structurally coherent region at a time
keeps the pooled distribution closer to unimodal and interpretable.

## Read the per-peptide diagnostic first

`summary` lists the matched per-peptide diagnostic first
(`n_peptides_tested`, `n_peptides_significant`,
`frac_peptides_significant`), then the pooled region test. The pooled
test is approximate: it treats the region's peptides as independent
samples although overlapping peptides share residues, and it compares
peptide means without their replicate variance. For peptides measured in
every state,
[`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md)
(a two-way ANOVA with peptide as a block) is the more powerful
region-level test.

## Pooled significance is not peptide significance

The pooled test above answers "does this region's *typical* exchange
shift" - it says nothing about whether *some* peptides in the region
shift while others don't. A real, spatially localized effect (e.g. a
binding interface covering a handful of peptides inside an otherwise
unaffected 200-residue domain) gets averaged against every unaffected
peptide riding along in the same region and can come back completely
non-significant even though individual peptides in it clearly differ. To
surface this, `summary` also reports `n_peptides_tested`,
`n_peptides_significant`, and `frac_peptides_significant`: for every
peptide identity (Start, End, Charge and Sequence - each charge state
counts as its own peptide) shared between that pair's two states in that
region, a per-peptide matched Welch t-test (replicates of that one
peptide, state1 vs. state2 - the same kind of test
[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
and
[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
use), with p-values adjusted across the peptides of each pair of states,
and called significant when `p_adj < pvalue` and the difference exceeds
the critical interval. A region with a flat/non-significant pooled
result but `frac_peptides_significant > 0` means real, localized signal
exists in there - shrink `regions` to isolate it, or use peptide-matched
testing directly instead of pooling.

## What gets compared

`value`: `filepath`'s raw `# Deut` column (Da uptake, read as `X..Deut`)
is used directly by default (`value = "raw"`) - no per-peptide size
correction of any kind. Two other options exist, off by default:

- `"raw"`:

  Default. Raw deuterium uptake in Da, exactly the `# Deut` column from
  the input file - the most direct reading of the data, and what
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)/[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)/[`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md)
  use as well (their `percent = FALSE` default).

- `"percent"`:

  HDX Examiner's own `Deut %` column - uptake already expressed relative
  to that peptide's theoretical maximum.

- `"frac_deut"`:

  This package's own derived metric: `# Deut` divided by the peptide's
  number of exchangeable protons (see
  [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)) -
  a finer-grained per-residue normalization than `"percent"`. Useful
  specifically when pooling peptides of very different lengths together
  in one region and the *length* difference itself would otherwise
  dominate the pooled comparison (a longer peptide mechanically carries
  a larger raw Da value regardless of its intrinsic exchange rate) -
  trades that off against normalizing away part of the signal, which is
  why it is opt-in, not the default.

The per-peptide p-values do not depend on this choice (dividing one
peptide's two states by the same constant does not change a t-test's
p-value), but the critical-interval half of the per-peptide rule and the
*pooled* result (`shift_hl`/`p_value`/`p_omnibus`, which pools values
across different peptides) are computed on the chosen scale.

## `method`

- `"peptide"` (default):

  Each peptide contributes one value per state (its replicates averaged
  first). Recommended default - avoids replicate pseudo-replication
  (replicates of one peptide are not independent draws of "how this
  state exchanges"). Omnibus: Kruskal-Wallis. Pairwise: Wilcoxon
  rank-sum, reported alongside the Hodges-Lehmann shift estimate.

- `"classic"`:

  Pools every replicate row as its own observation (`replicates` x
  n_peptides points per state) - reproduces the pooling this package's
  own prior ad hoc analyses used. Triggers a warning: replicates within
  a peptide are highly correlated, so this inflates the effective sample
  size and understates p-values. Kept for reproducing prior results, not
  recommended for new analyses.

- `"resample"`:

  Peptide-level like `"peptide"`, but p-values come from a
  label-permutation test (permuting which peptide belongs to which state
  and recomputing the test statistic) rather than the
  Kruskal-Wallis/Wilcoxon asymptotic approximation, and the pairwise
  `shift_hl` gets a percentile bootstrap confidence interval (resampling
  peptides with replacement within each state). Slower (`n_perm` +
  `n_boot` resamples per pair/region); more defensible with small
  peptide counts per state where the asymptotic tests are unreliable.

## See also

[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md),
[`all_data_long_input()`](https://mkajano.github.io/HDXBoxeR/reference/all_data_long_input.md),
[`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md),
[`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md),
[`plot_coverage_diagnostics()`](https://mkajano.github.io/HDXBoxeR/reference/plot_coverage_diagnostics.md)

## Examples

``` r
# \donttest{
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
res <- compare_states_unmatched(file_nm, times = "1800.00s",
                                 regions = list(nterm = c(1, 100), cterm = c(101, 225)))
res$summary
#>   region  state1 state2 n_peptides_tested n_peptides_significant
#> 1  nterm Unbound  bound                 9                      5
#> 2  cterm Unbound  bound                19                      6
#>   frac_peptides_significant n1 n2   shift_hl ci_low ci_high   p_value p_omnibus
#> 1                 0.5555556 11  9 -1.0886667     NA      NA 0.2874954 0.2706274
#> 2                 0.3157895 21 22 -0.3056667     NA      NA 0.5762709 0.5680073
#>       p_adj
#> 1 0.2874954
#> 2 0.5762709
res$coverage
#>   region   state n_peptides mean_length frac_shared_all_states
#> 1  nterm Unbound         11    15.72727              0.8181818
#> 2  nterm   bound          9    14.22222              0.8181818
#> 3  cterm Unbound         21    12.33333              0.8695652
#> 4  cterm   bound         22    12.09091              0.8695652
# }
```
