# Compare peptide-level HDX uptake across three or more protein states

Every other significance function in this package
([`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
[`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md),
[`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md))
compares exactly two states at a time. This function does two things for
three or more states: (1) an omnibus test per peptide across *all*
requested states (does this peptide differ across any of them at all) -
by default a Welch one-way ANOVA
([`stats::oneway.test()`](https://rdrr.io/r/stats/oneway.test.html),
unequal variances), reported but not used to gate significance (see
Details) - and (2) the package's usual Hageman & Weis hybrid test
(pooled
[`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md) +
per-pair Welch p-value) for every pair of states (or every `reference`
vs. the rest, if `reference` is given). A pair is significant when its
Welch p-value (raw by default, or adjusted with `p_adjust` across the
peptides of that pair) is below `pvalue` AND its difference exceeds the
pooled CI - the same rule as every other function in the package.

## Usage

``` r
compare_states_tp(
  filepath,
  times,
  states = NULL,
  by = c("sequence", "position"),
  reference = NULL,
  regions = NULL,
  replicates = NULL,
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  omnibus = c("welch", "kruskal"),
  percent = FALSE,
  p_adjust_method
)
```

## Arguments

- filepath:

  filepath to the input file (same format as
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)/[`nm_states()`](https://mkajano.github.io/HDXBoxeR/reference/nm_states.md)
  expect) - this function calls
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  internally, once per pair plus once for all states together, rather
  than taking an already-built wide table, precisely because that wide
  format has already discarded the availability information
  available-case matching needs (see Details).

- times:

  a single `Deut.Time` value to compare at.

- states:

  character vector of three or more states to compare. Default `NULL`
  uses every state found in `filepath`
  ([`nm_states()`](https://mkajano.github.io/HDXBoxeR/reference/nm_states.md)).

- by:

  matching rule for identifying "the same peptide" across states.
  `"sequence"` (default) requires identical Sequence, Start, End and
  Charge in both states
  ([`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)'s
  `seq_match = FALSE`). `"position"` matches on Start, End and Charge
  only, ignoring Sequence
  ([`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)'s
  `seq_match = TRUE`), so a peptide covering a point mutation is still
  compared; its `Sequence` is then reported as the per-state sequences
  joined by `"."`. In both modes `peptide_id` is `Start_End_zCharge`
  (plus the sequence in `"sequence"` mode), so each charge state is its
  own peptide.

- reference:

  optional state name. If given, only `reference` vs. each other state
  is tested (N-1 pairs, Dunnett-style), instead of every pair
  (N(N-1)/2) - both fewer simultaneous tests and usually the
  actually-relevant contrast for a wild-type/apo vs. several
  perturbations design. Default `NULL` tests every pair.

- regions:

  optional named list of `c(start, end)` residue ranges (e.g.
  `list(NTR = c(1, 85), ACD = c(86, 165))`). When supplied, only
  peptides fully contained in a region (`Start >= start & End <= end`)
  are compared, run separately per region, and the result gets a
  `region` column. Default `NULL` compares every peptide with no
  filtering (no `region` column in the output).

- replicates:

  optional maximum number of replicates per state (passed to
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)).
  Default `NULL` uses every replicate; states may differ.

- alpha:

  alpha for the pooled critical interval
  ([`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md)),
  applied per pair. Default `0.01`.

- pvalue:

  significance cutoff applied to the per-pair p-value (raw, or adjusted
  with `p_adjust`). Default `0.01`.

- p_adjust:

  multiple-testing correction passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw Welch p-values, as in every other function) or e.g.
  `"BH"`. Pairwise p-values are adjusted across the peptides of each
  pair of states (the same family as everywhere in the package); the
  omnibus p-values are adjusted across peptides.

- omnibus:

  omnibus test across all states per peptide: `"welch"` (default, Welch
  one-way ANOVA) or `"kruskal"` (Kruskal-Wallis).

- percent:

  passed to
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md) -
  `FALSE` (default) compares uptake in Daltons, `TRUE` percent
  deuteration, `"FD"` percent of the fully deuterated control,
  `"per_residue"` uptake per exchangeable amide (removes peptide length
  from the critical-interval condition).

- p_adjust_method:

  deprecated; earlier versions adjusted with `p_adjust_method = "BH"` by
  default. Supplying it sets `p_adjust`, with a warning.

## Value

data.frame with one row per peptide-by-state-pair (by region, if
`regions` is supplied): `peptide_id`, `Start`, `End`, `Sequence`,
`state1`, `state2`, `delta` (`state2 - state1`), `p_value` (raw per-pair
Welch p-value), `ci` (pooled CI for that pair), `p_adj` (per-pair
p-value after `p_adjust`; equal to `p_value` when `"none"`),
`p_omnibus`/`p_omnibus_adj` (`NA` for peptides absent from at least one
requested state), and `is_sig`.

## Which replicate SD feeds the pooled CI

Each state's own raw replicate SD, pooled via
[`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md) -
the same threshold as
[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)'s
default (`sd_source = "raw"`),
[`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md),
[`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md)
and the interactive plots, so `ci` values here match those functions for
the same pair (when those are given a table of the same two states; on a
table of 3+ states they pool over the peptides found in every state, so
the CI can differ slightly). It treats every state symmetrically and
needs no reference state.

## Available-case, not complete-case, matching

[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)'s
wide-format merge is an inner join across every state passed to it - a
peptide missing from *any one* of N states is dropped entirely, so
calling
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
once with all N states would silently shrink the comparable peptide set
combinatorially as N grows (a peptide identified in 90% of runs per
state survives in only 0.9^N of cases). This function instead calls
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
**separately per pair**, so a peptide only needs to be present in the
*two* states of that specific pair to be tested there - "available-case"
matching. The one place complete-case is unavoidable is the omnibus test
itself: it cannot test a peptide against a group it has no data for, so
the omnibus p-value is computed from the peptide set present in *all*
requested states (via one additional
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
call with every state), and peptides outside that set simply get `NA` in
`p_omnibus`/`p_omnibus_adj`.

## The omnibus test

The default omnibus test is a Welch one-way ANOVA, which keeps useful
power at the small replicate counts typical of HDX (3 per state). The
rank-based Kruskal-Wallis test (`omnibus = "kruskal"`) is available but
is very coarse at this size: with 3 states of 3 replicates its smallest
attainable p-value is about 0.0036. Either way the omnibus p-value is
reported for context (`p_omnibus`, and `p_omnibus_adj` adjusted across
peptides with `p_adjust`); `is_sig` depends only on the per-pair rule.

## See also

[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
[`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md),
[`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md),
[`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)

## Examples

``` r
# \donttest{
file_3state <- system.file("extdata", "three_states.csv", package = "HDXBoxeR")
res <- compare_states_tp(file_3state, times = "1800.00s", reference = "WT_pH7.5")
#> 22 peptide row(s) dropped: 22 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4 R120G_ph7.5=4
#> no csv file written
#> 10 peptide row(s) dropped: 10 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4
#> no csv file written
#> 14 peptide row(s) dropped: 14 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 R120G_ph7.5=4
#> no csv file written
head(res)
#>                        peptide_id Start End                Sequence   state1
#> 1                1_9_z2_MDIAIHHPW     1   9               MDIAIHHPW WT_pH7.5
#> 2           1_13_z3_MDIAIHHPWIRRP     1  13           MDIAIHHPWIRRP WT_pH7.5
#> 3 1_23_z3_MDIAIHHPWIRRPFFPFHSPSRL     1  23 MDIAIHHPWIRRPFFPFHSPSRL WT_pH7.5
#> 4 1_23_z4_MDIAIHHPWIRRPFFPFHSPSRL     1  23 MDIAIHHPWIRRPFFPFHSPSRL WT_pH7.5
#> 5 1_23_z5_MDIAIHHPWIRRPFFPFHSPSRL     1  23 MDIAIHHPWIRRPFFPFHSPSRL WT_pH7.5
#> 6  2_23_z4_DIAIHHPWIRRPFFPFHSPSRL     2  23  DIAIHHPWIRRPFFPFHSPSRL WT_pH7.5
#>        state2   delta      p_value        ci    p_omnibus        p_adj
#> 1 D109H_pH7.5 0.09850 7.793244e-03 0.1860082 1.063243e-02 7.793244e-03
#> 2 D109H_pH7.5 0.24650 9.561999e-04 0.1860082 1.420587e-03 9.561999e-04
#> 3 D109H_pH7.5 0.73650 2.157453e-05 0.1860082 4.837629e-05 2.157453e-05
#> 4 D109H_pH7.5 0.68750 9.520933e-04 0.1860082 1.254377e-03 9.520933e-04
#> 5 D109H_pH7.5 0.80325 1.477489e-04 0.1860082 2.465141e-05 1.477489e-04
#> 6 D109H_pH7.5 0.82500 5.675467e-04 0.1860082 1.812711e-04 5.675467e-04
#>   p_omnibus_adj is_sig
#> 1  1.063243e-02  FALSE
#> 2  1.420587e-03   TRUE
#> 3  4.837629e-05   TRUE
#> 4  1.254377e-03   TRUE
#> 5  2.465141e-05   TRUE
#> 6  1.812711e-04   TRUE
# }
```
