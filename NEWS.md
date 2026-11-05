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
