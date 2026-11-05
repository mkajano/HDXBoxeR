# Changelog

## HDXBoxeR 0.1.0

### New features

- [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md):
  state 1 vs state 2 deuteration plot with per-peptide Welch confidence
  ellipses (an ellipse touching the diagonal means the peptide’s own
  p-value is above the cutoff) and the pooled Hageman & Weis
  critical-interval band.
- Interactive, linked views:
  [`interactive_volcano()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_volcano.md),
  [`interactive_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_vent.md),
  [`interactive_woods()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_woods.md),
  [`interactive_fissure_explorer()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_explorer.md)
  (plotly/crosstalk).
- Multi-state comparisons:
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
  (three or more states, Welch ANOVA omnibus plus pairwise hybrid tests,
  optional BH correction; `by = "sequence"` or `by = "position"` peptide
  matching) and
  [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
  (region-level comparison when peptide maps differ between states, with
  a per-peptide diagnostic), with
  [`plot_compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_tp.md),
  [`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md)
  and
  [`plot_coverage_diagnostics()`](https://mkajano.github.io/HDXBoxeR/reference/plot_coverage_diagnostics.md).
- Per-residue analysis:
  [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md),
  [`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md),
  [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md),
  [`residue_regions_subset_summary()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions_subset_summary.md).
- Structure mapping:
  [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md),
  [`plot_hdx_structure_legend()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure_legend.md).
- Data quality:
  [`plot_peptide_coverage_map()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_coverage_map.md);
  long-format helpers
  [`boxer_input_to_long()`](https://mkajano.github.io/HDXBoxeR/reference/boxer_input_to_long.md)
  and
  [`all_data_long_input()`](https://mkajano.github.io/HDXBoxeR/reference/all_data_long_input.md).
- New vignettes: fissure vent, interactive explorer, multi-state
  (anova), structure mapping, data quality. The unfinished
  correlation-plots vignette was removed.

### Interface and documentation

- Argument names are consistent: `pvalue` (was `pv_cutoff`) and `col`
  (was `cola` / `colors_initial`). The old names still work, with a
  deprecation warning. As in 0.0.2, `pv_cutoff` sets both the p-value
  cutoff and the critical-interval `alpha` (unless `alpha` is given).
- The plotting and PyMOL functions now have separate `alpha` and
  `pvalue` arguments, placed where `pv_cutoff` was. Calls that pass
  these values by position
  (e.g. `plot_heat_map_tp(a, 3.5, ranges, 0.05, 3)`) now set different
  arguments: name the arguments. `CI_factor` was removed from the robot
  and Woods plot functions; use `alpha`.
- PyMOL script functions gain `resno_offset` to match structure
  numbering.
- Internal helper functions are marked internal (still exported; hidden
  from the help index).
- The `_proc` functions are documented as percent-deuteration
  counterparts of their uptake versions (significance on uptake,
  differences shown in percent deuteration).
- New example data `three_states.csv` (wild type and two point mutants,
  uneven replicates; Woods et al. 2025, PNAS).
- New
  [`vignette("input_formats")`](https://mkajano.github.io/HDXBoxeR/articles/input_formats.md);
  README rewritten; pkgdown reference grouped by task; main vignette
  shows its summary/output code and points to the newer tools; stale
  references to removed functions fixed.

### Reporting

- New
  [`render_report()`](https://mkajano.github.io/HDXBoxeR/reference/render_report.md):
  one-command HTML report with the analysis parameters, the summary
  table, data-quality checks, significant-peptide counts and
  fissure-vent plots per time point, and optionally region ANOVAs and
  whole-time-course results; the tidy results table is written next to
  it.
- New
  [`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md):
  the summary table recommended by Masson et al.
  2019. - computed items (time course, controls, replicates, peptides,
          coverage, length, redundancy, back-exchange, repeatability,
          significance threshold) plus your own experimental details via
          `metadata`; csv or markdown.
- New
  [`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md):
  tidy per-peptide, per-time, per-pair results (means, SDs, n,
  difference, p-value, adjusted p-value, CI, call), written as csv with
  a `_parameters.txt` file holding the parameters and
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
- Results of
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md),
  [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md),
  [`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md),
  [`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md),
  [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md),
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  and the new tables carry their parameters, the package version and the
  date as attribute `"hdxboxer"`.
- [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
  gains `seed` for reproducible resampling.
- [`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md):
  `Mean.Peptide_Mass_[Da]` is now the measured neutral peptide mass
  (Charge x (centroid m/z - proton mass)); it previously held the
  centroid m/z. The percent column is labelled `[%]` instead of `[Da]`,
  column names are no longer mangled (`_[Da]` instead of `_.Da.`), and
  the table is returned invisibly.
- [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md)/[`all_summary()`](https://mkajano.github.io/HDXBoxeR/reference/all_summary.md):
  average peptide length counts residues (End - Start + 1); it was one
  residue too short.

### Statistics: one significance rule, regions and time courses

- One significance rule everywhere: a peptide differs when its Welch
  p-value is below `pvalue` AND its difference exceeds the pooled
  critical interval. p-values are raw by default; every function that
  calls significance now takes `p_adjust` (e.g. `"BH"`), including
  [`pv_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timepoint.md),
  [`pv_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timecourse.md),
  the heat map, volcano, Woods, robot and peptide plots, the PyMOL
  scripts,
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
  the interactive plots,
  [`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md)
  and
  [`verbose_timepoint_output()`](https://mkajano.github.io/HDXBoxeR/reference/verbose_timepoint_output.md).
- [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md):
  raw p-values by default (`p_adjust = "none"`; previously BH) and a
  Welch one-way ANOVA as the reported omnibus test
  (`omnibus = "kruskal"` for the previous Kruskal-Wallis test).
  `p_adjust_method` is deprecated. A peptide in two overlapping regions
  now enters the multiple-testing correction once.
- [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md):
  peptides are called by the package-wide rule (previously BH-adjusted p
  \< `fdr_cutoff` with no critical interval; `fdr_cutoff` still works,
  deprecated, and reproduces that rule). Charge states of one peptide
  count once towards `min_peptides`, and the per-residue effect is the
  difference per exchangeable amide.
- [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md):
  the per-peptide diagnostic uses the package-wide rule, is listed first
  in `summary`, and its per-peptide results are returned as `peptides`;
  the pooled test is documented as approximate. `p_adjust` replaces
  `p_adjust_method` (default now `"none"`).
- New
  [`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md):
  region-level two-way ANOVA (peptide + state) for matched peptides,
  with each peptide as its own control; pairwise state contrasts, a
  state x peptide heterogeneity test, and an optional model over all
  time points with a state x time test.
- New
  [`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md):
  per-peptide summed difference over the time course with a global
  critical interval, and the number of time points at which a peptide is
  significant.
  [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md),
  [`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
  and
  [`robot_plot_All()`](https://mkajano.github.io/HDXBoxeR/reference/robot_plot_All.md)
  gain `min_timepoints`.
- New `percent = "per_residue"` in
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)/[`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md):
  uptake per exchangeable amide, so peptide length does not affect the
  critical-interval half of the rule.
- [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  returns (invisibly) the table of what it drew.
- The multi-state vignette is now
  [`vignette("multistate")`](https://mkajano.github.io/HDXBoxeR/articles/multistate.md).

### Fully deuterated control and quality checks

- [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)/[`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md)
  gain `percent = "FD"`: percent of the measured fully deuterated
  control, `100 * uptake / FD uptake` of the same peptide and state (FD
  replicates averaged, charge states pooled), which corrects for
  back-exchange and D2O content. Peptides without a usable FD value
  (missing, or FD uptake \<= 0) keep their row with empty (`NA`) values.
  `percent = TRUE` still returns the input’s own percent deuteration
  (HDXExaminer `Deut %`, relative to the theoretical maximum).
- New
  [`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md)
  lists per peptide and state: missing undeuterated control, uptake
  reported as 0 at every time (a missing reference rather than a real
  measurement), negative uptake, FD uptake, theoretical maximum, FD
  fraction, and FD missing / \<= 0 / above the maximum. It only reports;
  nothing is removed.
- [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)
  computes default percent deuteration against the same theoretical
  maximum as HDXExaminer and DynamX (length - 1 - prolines after the
  first residue).

### Input formats

- New
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  reads HDXExaminer “All results” tables and Waters DynamX cluster and
  state exports (DynamX 2.x and 3.x), detecting the format from the
  header.
  [`read_hdxexaminer()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdxexaminer.md)
  and
  [`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md)
  read one format explicitly. For DynamX cluster data each raw file is a
  replicate and uptake is computed from the centroid mass; on a public
  matched cluster/state data set it agrees with DynamX’s own uptake (r =
  0.9999). `fd =` marks the fully deuterated control (a control state,
  or each state’s last time point); `combine_charges = TRUE` merges
  charge states.
- New
  [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)
  builds HDXBoxeR input from a long table exported by any other
  software, with times in s, min or h.
- Every function that takes `filepath` now also accepts a DynamX file or
  a table from
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)/[`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md);
  reading once and passing the table avoids re-reading the file in every
  call.

### Uneven and missing replicates

- [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  and
  [`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md)
  now use every replicate of every state by default, so states may have
  different numbers of replicates (e.g. 3 and 4). Previously all states
  were cut to the smallest count and the extra replicates were silently
  dropped; `replicates = n` still caps the number.
- A peptide with a missing measurement is kept with `NA` instead of
  being dropped, as long as it has at least `min_replicates` (new
  argument, default 2) values in every state. Dropped peptides are
  reported.
- [`sd_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timepoint.md),
  [`ave_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/ave_timepoint.md),
  [`pv_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timepoint.md),
  [`sd_ci_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_ci_timepoint.md)
  and
  [`pv_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timecourse.md)
  find each state’s replicate columns by name and ignore missing values;
  [`pv_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timepoint.md)
  returns `NA` instead of failing when a peptide has fewer than two
  values or identical replicates.
- [`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md)
  (and so
  [`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md),
  [`CI_tc()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tc.md) and
  every plot using them) pools each state’s SDs weighted by their
  degrees of freedom and uses each state’s own number of replicates.
  With equal, complete replicates the result is unchanged.
- The results of
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)/[`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md)
  record their replicates, states, times and units as attributes.
- State names are matched exactly before abbreviations are tried, so a
  state “WT” is no longer confused with “WT_mut”.
- `output_tc(seq_match = )` now means the same as in
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  (`FALSE`: identical sequence required; `TRUE`: match by position). The
  two functions previously used opposite meanings. With
  `seq_match = TRUE`,
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  keeps its usual column order (Sequence fourth).
- [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
  uses every replicate by default (`replicates = NULL`).

### Changes to results

- Robot plots, Woods plots,
  [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)
  and
  [`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md):
  p-values and averages came from two tables sorted differently (by
  Start as text in one, as a number in the other) and were combined row
  by row, so the wrong peptides could be called significant. This was
  also the case in 0.0.2. Both tables now have the same order.

- [`plots_vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_vol_tp.md)
  colors each state with its own critical interval instead of the
  largest one.

- [`plot_heat_map_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_tp_proc.md):
  the default `ranges` is `c(-Inf, seq(-30, 30, by = 10), Inf)`, as
  documented (it was ten times smaller).

- [`backHX_calculations()`](https://mkajano.github.io/HDXBoxeR/reference/backHX_calculations.md)
  averages the replicate centroids of each peptide before the
  calculation, so the range is over peptides; a state without FD or
  undeuterated data gives `NA`.

- [`output_UD()`](https://mkajano.github.io/HDXBoxeR/reference/output_UD.md)
  and the time-course tables treat a missing undeuterated uptake
  (HDXExaminer `n/a`, DynamX) as 0, so
  [`average_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/average_timecourse.md),
  [`sd_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timecourse.md)
  and their `_proc` versions work with such data.

- [`backHX_calculations()`](https://mkajano.github.io/HDXBoxeR/reference/backHX_calculations.md)
  now computes back-exchange from the measured centroid times charge,
  `(Exp Cent[FD] - Exp Cent[0s]) * Charge`. It previously used the
  theoretical centroid m/z without the charge factor, which overstated
  back-exchange roughly in proportion to the charge state.

- [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  (default `sd_source = "raw"`), the interactive plots and
  [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)
  now build the critical interval from each state’s raw replicate SD, so
  it matches
  [`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md),
  [`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md)
  and the heat map, Woods and robot plots. `sd_source = "CI"` keeps the
  previous, wider band.

- All critical-interval functions
  ([`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md),
  [`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md),
  [`CI_tc()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tc.md),
  [`CI_single()`](https://mkajano.github.io/HDXBoxeR/reference/CI_single.md),
  [`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md))
  use the same two-tailed t-value, `qt(alpha/2, df)`, as in HDXBoxeR
  0.0.2 and Janowska et al. (2024).

- Multiple-testing correction (`p_adjust`) uses one family everywhere:
  the peptides of one pair of states at one time point.

- [`sd_ci_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_ci_timepoint.md)
  uses the correct conversion from a two-sided interval to an SD (factor
  of 2 and `alpha/2`).

- [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md)
  reports the real number of replicates.

### Bug fixes

- `replicates` is now honoured by
  [`woods_CI_plot()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot.md),
  [`woods_CI_plot_frac()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot_frac.md),
  [`robot_plot_All()`](https://mkajano.github.io/HDXBoxeR/reference/robot_plot_All.md),
  [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md),
  [`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md),
  [`pymol_script_significant_residue()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_residue.md),
  [`pymol_script_significant_residue_proc()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_residue_proc.md)
  and
  [`summary_sd_CI()`](https://mkajano.github.io/HDXBoxeR/reference/summary_sd_CI.md).
- [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)/[`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
  keep significant peptides from every state, not only the last one.
- [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
  treats each charge state as its own peptide instead of pooling charge
  states as extra replicates.
- [`plot_cor_pts_elps()`](https://mkajano.github.io/HDXBoxeR/reference/plot_cor_pts_elps.md):
  ellipse position and axis fixed.
- [`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md):
  fully deuterated rows no longer report a fake SD of 100.
- [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)/[`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md)
  use the same exchangeable-residue convention as
  [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md).
- [`all_data_long_input()`](https://mkajano.github.io/HDXBoxeR/reference/all_data_long_input.md)
  no longer fails on an undefined variable;
  [`boxer_input_to_long()`](https://mkajano.github.io/HDXBoxeR/reference/boxer_input_to_long.md)
  parses percent-deuteration columns correctly.
- [`vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/vol_tp.md)
  default palette size no longer reads
  [`stats::df`](https://rdrr.io/r/stats/Fdist.html).
- Peptides are matched between states without the retention time, so a
  small RT shift no longer splits one peptide into two rows.
- [`output_FD()`](https://mkajano.github.io/HDXBoxeR/reference/output_FD.md),
  [`output_UD()`](https://mkajano.github.io/HDXBoxeR/reference/output_UD.md)
  and their `_proc` versions return every state (previously only the
  last one);
  [`output_UD()`](https://mkajano.github.io/HDXBoxeR/reference/output_UD.md)
  also finds an undeuterated control written as `"0s"`.
- [`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
  returns the peptide keys (Start, End, Sequence, Charge);
  [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)
  returns rows of the table used by the robot plots.
- [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md)
  and
  [`all_summary()`](https://mkajano.github.io/HDXBoxeR/reference/all_summary.md):
  coverage is computed over `protein_range` (new argument), and
  redundancy over covered residues.
  [`all_summary()`](https://mkajano.github.io/HDXBoxeR/reference/all_summary.md)
  works without back-exchange data.
- [`summary_sd_CI()`](https://mkajano.github.io/HDXBoxeR/reference/summary_sd_CI.md)
  pools SDs correctly with uneven replicates.
- The `_proc` functions stop when the uptake and percent tables are not
  aligned row by row.
- PyMOL scripts take the state names from the column names correctly.
- [`significant_peptide_uptake()`](https://mkajano.github.io/HDXBoxeR/reference/significant_peptide_uptake.md)
  no longer fails on missing values.
- [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  gains `xlim`/`ylim`; its returned table includes `Charge`.
- `read_dynamx(fd = "last")` warns when states have different longest
  exposures.
- [`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md):
  a state without its own undeuterated run (e.g. an FD control) uses the
  undeuterated mass of the same peptide from the other states, with a
  message.
- [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
  stops when given more than one time point.
- [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
  and
  [`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md)
  use the plain sequence and the Start-End range, so joined sequences
  (`seq_match = TRUE`) or DynamX modifications no longer spread uptake
  beyond the peptide;
  [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
  ignores modifications.
- PyMOL functions work with the default `path = ""` (working directory).
- [`pymol_script_significant_residue_proc()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_residue_proc.md)
  keeps `input_up` first (as in 0.0.2) and swaps the two tables, with a
  warning, when they are given the other way round.
- Interactive plots treat charge states of one peptide as separate
  peptides; vector `alpha`/`pvalue` no longer break
  [`interactive_woods()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_woods.md)
  and
  [`significant_peptides_table()`](https://mkajano.github.io/HDXBoxeR/reference/significant_peptides_table.md).
- [`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md)
  uses every replicate by default, handles uneven replicates and missing
  centroids.
- [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)
  matches state names exactly.
- State names in plot legends and titles are no longer cut short for
  percent tables.
- Average, difference and volcano plots work with `NA` values (e.g.
  `percent = "FD"` tables).
- `seq_match = TRUE` with three or more states lists each sequence once;
  time labels in minutes or hours sort by time.
- [`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md)
  gains `states`; the report summary uses the report’s states.

### Documentation and packaging

- Test suite added (testthat 3rd edition).
- [`?All_results_table`](https://mkajano.github.io/HDXBoxeR/reference/All_results_table.md)
  documents the example input file.
- `URL` and `BugReports` added to DESCRIPTION.

## HDXBoxeR 0.0.2

CRAN release: 2024-08-24

- Version published on CRAN; see Janowska et al. (2024) Bioinformatics,
  <doi:10.1093/bioinformatics/btae479>.
