# HDXBoxeR (development version)

## New features

* `plot_fissure_vent()`: state 1 vs state 2 deuteration plot with per-peptide
  Welch confidence ellipses (an ellipse touching the diagonal means the
  peptide's own p-value is above the cutoff) and the pooled Hageman & Weis
  critical-interval band.
* Interactive, linked views: `interactive_volcano()`, `interactive_fissure_vent()`,
  `interactive_woods()`, `interactive_fissure_explorer()` (plotly/crosstalk).
* Multi-state comparisons: `compare_states_tp()` (three or more states,
  Kruskal-Wallis omnibus plus pairwise hybrid tests with BH correction;
  `by = "sequence"` or `by = "position"` peptide matching) and
  `compare_states_unmatched()` (region-level comparison when peptide maps
  differ between states, with a per-peptide diagnostic), with
  `plot_compare_states_tp()`, `plot_compare_states_unmatched()` and
  `plot_coverage_diagnostics()`.
* Per-residue analysis: `residue_deconvolution()`, `residue_level_average()`,
  `residue_regions()`, `residue_regions_subset_summary()`.
* Structure mapping: `plot_hdx_structure()`, `plot_hdx_structure_legend()`.
* Data quality: `plot_peptide_coverage_map()`; long-format helpers
  `boxer_input_to_long()` and `all_data_long_input()`.
* New vignettes: fissure vent, interactive explorer, multi-state (anova),
  structure mapping, data quality.

## Interface and documentation

* Argument names are consistent: `pvalue` (was `pv_cutoff`) and `col` (was
  `cola` / `colors_initial`). The old names still work, with a deprecation
  warning; positional calls are unaffected.
* PyMOL script functions gain `resno_offset` to match structure numbering.
* Internal helper functions are marked internal (still exported; hidden
  from the help index).
* The `_proc` functions are documented as percent-deuteration counterparts
  of their uptake versions (significance on uptake, differences shown in
  percent deuteration).
* New example data `three_states.csv` (wild type and two point mutants,
  uneven replicates; Woods et al. 2025, PNAS).
* New `vignette("input_formats")`; README rewritten; pkgdown reference
  grouped by task; main vignette shows its summary/output code and points
  to the newer tools; stale references to removed functions fixed.

## Reporting

* New `render_report()`: one-command HTML report with the analysis
  parameters, the summary table, data-quality checks, significant-peptide
  counts and fissure-vent plots per time point, and optionally region
  ANOVAs and whole-time-course results; the tidy results table is written
  next to it.
* New `hdx_summary_table()`: the summary table recommended by Masson et al.
  (2019) - computed items (time course, controls, replicates, peptides,
  coverage, length, redundancy, back-exchange, repeatability, significance
  threshold) plus your own experimental details via `metadata`; csv or
  markdown.
* New `hdx_results_table()`: tidy per-peptide, per-time, per-pair results
  (means, SDs, n, difference, p-value, adjusted p-value, CI, call), written
  as csv with a `_parameters.txt` file holding the parameters and
  `sessionInfo()`.
* Results of `compare_states_tp()`, `compare_states_unmatched()`,
  `compare_regions()`, `compare_timecourse()`, `residue_regions()`,
  `plot_fissure_vent()` and the new tables carry their parameters, the
  package version and the date as attribute `"hdxboxer"`.
* `compare_states_unmatched()` gains `seed` for reproducible resampling.
* `output_prep()`: `Mean.Peptide_Mass_[Da]` is now the measured neutral
  peptide mass (Charge x (centroid m/z - proton mass)); it previously held
  the centroid m/z. The percent column is labelled `[%]` instead of `[Da]`,
  column names are no longer mangled (`_[Da]` instead of `_.Da.`), and the
  table is returned invisibly.
* `general_info()`/`all_summary()`: average peptide length counts residues
  (End - Start + 1); it was one residue too short.

## Statistics: one significance rule, regions and time courses

* One significance rule everywhere: a peptide differs when its Welch
  p-value is below `pvalue` AND its difference exceeds the pooled
  critical interval. p-values are raw by default; every function that
  calls significance now takes `p_adjust` (e.g. `"BH"`), including
  `pv_timepoint()`, `pv_timecourse()`, the heat map, volcano, Woods, robot
  and peptide plots, the PyMOL scripts, `plot_fissure_vent()`, the
  interactive plots, `diff_peptides_prep()` and `verbose_timepoint_output()`.
* `compare_states_tp()`: raw p-values by default (`p_adjust = "none"`;
  previously BH) and a Welch one-way ANOVA as the reported omnibus test
  (`omnibus = "kruskal"` for the previous Kruskal-Wallis test).
  `p_adjust_method` is deprecated. A peptide in two overlapping regions now
  enters the multiple-testing correction once.
* `residue_regions()`: peptides are called by the package-wide rule
  (previously BH-adjusted p < `fdr_cutoff` with no critical interval;
  `fdr_cutoff` still works, deprecated, and reproduces that rule). Charge
  states of one peptide count once towards `min_peptides`, and the
  per-residue effect is the difference per exchangeable amide.
* `compare_states_unmatched()`: the per-peptide diagnostic uses the
  package-wide rule, is listed first in `summary`, and its per-peptide
  results are returned as `peptides`; the pooled test is documented as
  approximate. `p_adjust` replaces `p_adjust_method` (default now `"none"`).
* New `compare_regions()`: region-level two-way ANOVA (peptide + state) for
  matched peptides, with each peptide as its own control; pairwise state
  contrasts, a state x peptide heterogeneity test, and an optional model
  over all time points with a state x time test.
* New `compare_timecourse()`: per-peptide summed difference over the time
  course with a global critical interval, and the number of time points at
  which a peptide is significant. `robot_indexes()`, `robot_indexes_df()`
  and `robot_plot_All()` gain `min_timepoints`.
* New `percent = "per_residue"` in `output_tp()`/`output_tc()`: uptake per
  exchangeable amide, so peptide length does not affect the
  critical-interval half of the rule.
* `plot_fissure_vent()` returns (invisibly) the table of what it drew.
* The multi-state vignette is now `vignette("multistate")`.

## Fully deuterated control and quality checks

* `output_tp()`/`output_tc()` gain `percent = "FD"`: percent of the measured
  fully deuterated control, `100 * uptake / FD uptake` of the same peptide
  and state (FD replicates averaged, charge states pooled), which corrects
  for back-exchange and D2O content. Peptides without a usable FD value
  (missing, or FD uptake <= 0) keep their row with empty (`NA`) values.
  `percent = TRUE` still returns the input's own percent deuteration
  (HDXExaminer `Deut %`, relative to the theoretical maximum).
* New `peptide_qc()` lists per peptide and state: missing undeuterated
  control, uptake reported as 0 at every time (a missing reference rather
  than a real measurement), negative uptake, FD uptake, theoretical
  maximum, FD fraction, and FD missing / <= 0 / above the maximum. It only
  reports; nothing is removed.
* `as_hdx_long()` computes default percent deuteration against the same
  theoretical maximum as HDXExaminer and DynamX (length - 1 - prolines
  after the first residue).

## Input formats

* New `read_hdx()` reads HDXExaminer "All results" tables and Waters
  DynamX cluster and state exports (DynamX 2.x and 3.x), detecting the
  format from the header. `read_hdxexaminer()` and `read_dynamx()` read one
  format explicitly. For DynamX cluster data each raw file is a replicate
  and uptake is computed from the centroid mass; on a public matched
  cluster/state data set it agrees with DynamX's own uptake (r = 0.9999).
  `fd =` marks the fully deuterated control (a control state, or each
  state's last time point); `combine_charges = TRUE` merges charge states.
* New `as_hdx_long()` builds HDXBoxeR input from a long table exported by
  any other software, with times in s, min or h.
* Every function that takes `filepath` now also accepts a DynamX file or a
  table from `read_hdx()`/`as_hdx_long()`; reading once and passing the
  table avoids re-reading the file in every call.

## Uneven and missing replicates

* `output_tp()` and `output_tc()` now use every replicate of every state by
  default, so states may have different numbers of replicates (e.g. 3 and
  4). Previously all states were cut to the smallest count and the extra
  replicates were silently dropped; `replicates = n` still caps the number.
* A peptide with a missing measurement is kept with `NA` instead of being
  dropped, as long as it has at least `min_replicates` (new argument,
  default 2) values in every state. Dropped peptides are reported.
* `sd_timepoint()`, `ave_timepoint()`, `pv_timepoint()`, `sd_ci_timepoint()`
  and `pv_timecourse()` find each state's replicate columns by name and
  ignore missing values; `pv_timepoint()` returns `NA` instead of failing
  when a peptide has fewer than two values or identical replicates.
* `CI_2pts()` (and so `CI_tp()`, `CI_tc()` and every plot using them) pools
  each state's SDs weighted by their degrees of freedom and uses each
  state's own number of replicates. With equal, complete replicates the
  result is unchanged.
* The results of `output_tp()`/`output_tc()` record their replicates,
  states, times and units as attributes.
* State names are matched exactly before abbreviations are tried, so a
  state "WT" is no longer confused with "WT_mut".
* `output_tc(seq_match = )` now means the same as in `output_tp()`
  (`FALSE`: identical sequence required; `TRUE`: match by position). The
  two functions previously used opposite meanings. With `seq_match = TRUE`,
  `output_tp()` keeps its usual column order (Sequence fourth).
* `compare_states_tp()` uses every replicate by default (`replicates = NULL`).

## Changes to results

* `backHX_calculations()` now computes back-exchange from the measured
  centroid times charge, `(Exp Cent[FD] - Exp Cent[0s]) * Charge`. It
  previously used the theoretical centroid m/z without the charge factor,
  which overstated back-exchange roughly in proportion to the charge state.
* `plot_fissure_vent()` (default `sd_source = "raw"`), the interactive plots
  and `residue_regions()` now build the critical interval from each state's
  raw replicate SD, so it matches `CI_tp()`, `diff_peptides_prep()` and the
  heat map, Woods and robot plots. `sd_source = "CI"` keeps the previous,
  wider band.
* All critical-interval functions (`CI_2pts()`, `CI_tp()`, `CI_tc()`,
  `CI_single()`) use the same t-value, `qt(alpha, df)`.
* `sd_ci_timepoint()` uses the correct conversion from a two-sided interval
  to an SD (factor of 2 and `alpha/2`).
* `general_info()` reports the real number of replicates.

## Bug fixes

* `replicates` is now honoured by `woods_CI_plot()`, `woods_CI_plot_frac()`,
  `robot_plot_All()`, `robot_indexes()`, `robot_indexes_df()`,
  `pymol_script_significant_residue()`,
  `pymol_script_significant_residue_proc()` and `summary_sd_CI()`.
* `robot_indexes()`/`robot_indexes_df()` keep significant peptides from every
  state, not only the last one.
* `compare_states_unmatched()` treats each charge state as its own peptide
  instead of pooling charge states as extra replicates.
* `plot_cor_pts_elps()`: ellipse position and axis fixed.
* `output_prep()`: fully deuterated rows no longer report a fake SD of 100.
* `residue_deconvolution()`/`residue_level_average()` use the same
  exchangeable-residue convention as `nb_exch_deut()`.
* `all_data_long_input()` no longer fails on an undefined variable;
  `boxer_input_to_long()` parses percent-deuteration columns correctly.
* `vol_tp()` default palette size no longer reads `stats::df`.

## Documentation and packaging

* Test suite added (testthat 3rd edition).
* `?All_results_table` documents the example input file.
* `URL` and `BugReports` added to DESCRIPTION.

# HDXBoxeR 0.1.4

* Version published on CRAN; see Janowska et al. (2024) Bioinformatics,
  doi:10.1093/bioinformatics/btae479.
