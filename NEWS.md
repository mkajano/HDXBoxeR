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
