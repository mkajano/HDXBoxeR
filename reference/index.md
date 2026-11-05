# Package index

## Reading data

- [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  : Read HDX-MS data from HDXExaminer or DynamX
- [`read_hdxexaminer()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdxexaminer.md)
  : Read an HDXExaminer "All results" table
- [`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md)
  : Read a DynamX cluster or state data export
- [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)
  : Build HDXBoxeR input from a long table
- [`nm_states()`](https://mkajano.github.io/HDXBoxeR/reference/nm_states.md)
  : Lists names of states in data sets
- [`All_results_table`](https://mkajano.github.io/HDXBoxeR/reference/All_results_table.md)
  : Example input: All_results_table.csv from HDXExaminer
- [`three_states`](https://mkajano.github.io/HDXBoxeR/reference/three_states.md)
  : Example input: three protein states with point mutations

## Organizing data

Wide tables of replicate values used by every analysis.

- [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  : Prepares output for HDX-MS for the deuteration uptake or percent
  deuteration for the time points.
- [`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md)
  : Prepares output for HDX-MS for the deuteration uptake or percent
  deuteration for the time courses.
- [`select_indices()`](https://mkajano.github.io/HDXBoxeR/reference/select_indices.md)
  : Allows for selecting some peptide from input data
- [`boxer_input_to_long()`](https://mkajano.github.io/HDXBoxeR/reference/boxer_input_to_long.md)
  : Prepares the long input of the input data
- [`all_data_long_input()`](https://mkajano.github.io/HDXBoxeR/reference/all_data_long_input.md)
  : Raw all_data input
- [`output_FD()`](https://mkajano.github.io/HDXBoxeR/reference/output_FD.md)
  : Prepares output for HDX-MS Full deuteration data
- [`output_FD_proc()`](https://mkajano.github.io/HDXBoxeR/reference/output_FD_proc.md)
  : Prepares output for HDX-MS Full deuteration data for procent
  deuteration.
- [`output_UD()`](https://mkajano.github.io/HDXBoxeR/reference/output_UD.md)
  : Prepares output for HDX-MS Undeuterated sample data.
- [`output_UD_proc()`](https://mkajano.github.io/HDXBoxeR/reference/output_UD_proc.md)
  : Prepares output for HDX-MS Undeuterated data for procent
  deuteration.

## Data quality and summaries

- [`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md)
  : Per-peptide quality checks of the input data
- [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md)
  : Provides summary table for all data.sets.
- [`backHX_calculations()`](https://mkajano.github.io/HDXBoxeR/reference/backHX_calculations.md)
  : Summary of backexchange summary
- [`summary_sd_CI()`](https://mkajano.github.io/HDXBoxeR/reference/summary_sd_CI.md)
  : Provides summary table with Critical interval and standard deviation
  within the set.
- [`all_summary()`](https://mkajano.github.io/HDXBoxeR/reference/all_summary.md)
  : Returns full summary table.
- [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
  : Number of exchangeable protons
- [`plot_peptide_coverage_map()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_coverage_map.md)
  : Plot peptide coverage map with per-residue redundancy

## Statistics

Averages, SDs, Welch p-values and the Hageman & Weis critical interval.

- [`ave_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/ave_timepoint.md)
  : Returns average value for either uptake of procent data.
- [`sd_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timepoint.md)
  : Returns standard deviation for dataframe.
- [`pv_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timepoint.md)
  : Calculation of pvalue between first protein state and any other
  state from all_states file
- [`dif_ave()`](https://mkajano.github.io/HDXBoxeR/reference/dif_ave.md)
  : Returns data frame with difference of averages between State1 and
  other states provided.
- [`sd_ci_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_ci_timepoint.md)
  : Calculate what would be the standard deviation based on the CI
- [`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md)
  : Global confidence interval treshold from experimental standard
  deviation for 2 samples.
- [`CI_tp()`](https://mkajano.github.io/HDXBoxeR/reference/CI_tp.md) :
  Global confidence interval treshold from experimental standard
  deviation
- [`CI_single()`](https://mkajano.github.io/HDXBoxeR/reference/CI_single.md)
  : Global confidence interval treshold from experimental standard
  deviation for 1 sample
- [`diff_peptides_prep()`](https://mkajano.github.io/HDXBoxeR/reference/diff_peptides_prep.md)
  : Peptides significantly different between states

## Comparing states

- [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  : Fissure_vent plot

- [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
  : Compare peptide-level HDX uptake across three or more protein states

- [`plot_compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_tp.md)
  :

  Plot per-peptide deltas from
  [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md),
  by state pair and region

- [`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md)
  : Region-level ANOVA: compare protein states within residue regions

- [`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md)
  : Compare two states over a whole time course

- [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)
  : Compare pooled peptide-level HDX uptake across states without
  matching peptides

- [`plot_compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/plot_compare_states_unmatched.md)
  :

  Plot the pooled distributions behind
  [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)

- [`plot_coverage_diagnostics()`](https://mkajano.github.io/HDXBoxeR/reference/plot_coverage_diagnostics.md)
  :

  Plot the peptide-map coverage diagnostics from
  [`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)

- [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)
  : Cluster significant peptides into contiguous destabilized/protected
  regions

- [`residue_regions_subset_summary()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions_subset_summary.md)
  : Roll up derived regions against named residue subsets

## Per-residue values

- [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
  : Per-residue deconvolution of HDX-MS uptake from overlapping peptides
- [`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md)
  : Per-residue uptake as a length-weighted mean across covering
  peptides

## Plots - time points

Uptake (`*_tp`) and percent-deuteration (`*_proc`) versions; the `_proc`
functions show differences in percent deuteration and test significance
on the uptake table.

- [`plots_av_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_av_tp.md)
  : Returns average deuteration plot for timepoints in the data frame
- [`plots_av_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plots_av_tp_proc.md)
  : Returns average procent deuteration plot for time points
- [`plots_diff_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_diff_tp.md)
  : Returns difference in average plot for timepoints in the data frame
- [`plots_diff_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plots_diff_tp_proc.md)
  : Returns difference in average procent deuteration plot for
  timepoints in the data frame
- [`plots_vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_vol_tp.md)
  : Returns volcano plots for timepoints in the data frame
- [`plot_heat_map_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_tp.md)
  : Plots heat maps for significant peptides.
- [`plot_heat_map_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_tp_proc.md)
  : Plots heat maps for significant peptides.
- [`plot_heat_map_max_uptake_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_max_uptake_tp.md)
  : Plots heat maps for maximum uptake per residue.
- [`plot_heat_map_max_uptake_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_max_uptake_tp_proc.md)
  : Plots heat maps for maximum procent deuteration per residue.
- [`plot_peptide_sig_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_sig_tp.md)
  : Significant peptide plots.
- [`plot_peptide_sig_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_sig_tp_proc.md)
  : Draws peptides with significant difefrences between sets.
- [`boxplot_tp()`](https://mkajano.github.io/HDXBoxeR/reference/boxplot_tp.md)
  : Plots boxplots for all the averages in the set
- [`deuteration_woods_timepoints()`](https://mkajano.github.io/HDXBoxeR/reference/deuteration_woods_timepoints.md)
  : Return woods plots for the timepoints
- [`deuteration_woods_timepoints_frac()`](https://mkajano.github.io/HDXBoxeR/reference/deuteration_woods_timepoints_frac.md)
  : Return woods plots for the timepoints
- [`deuteration_woods_timepoints_frac_list()`](https://mkajano.github.io/HDXBoxeR/reference/deuteration_woods_timepoints_frac_list.md)
  : Return woods plots for the timepoints
- [`plot_cor_pts_elps()`](https://mkajano.github.io/HDXBoxeR/reference/plot_cor_pts_elps.md)
  : Correlation functions plots and prep data
- [`plot_correlation_functions()`](https://mkajano.github.io/HDXBoxeR/reference/plot_correlation_functions.md)
  : Plot correlated states

## Plots - time courses

- [`plots_av_tcourse()`](https://mkajano.github.io/HDXBoxeR/reference/plots_av_tcourse.md)
  : Generates average deuteration plot for the time-course.
- [`plot_heat_map_tc()`](https://mkajano.github.io/HDXBoxeR/reference/plot_heat_map_tc.md)
  : Plots heat maps for time courses.
- [`uptake_plots()`](https://mkajano.github.io/HDXBoxeR/reference/uptake_plots.md)
  : Uptake plots
- [`robot_plot_All()`](https://mkajano.github.io/HDXBoxeR/reference/robot_plot_All.md)
  : Returns a robot plot for comparisons of the timepoints samples
- [`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)
  : Returns indexes for peptides with significant difference between two
  sets
- [`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
  : Returns dataframe with peptides which exhibit significant difference
  between two sets
- [`woods_CI_plot()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot.md)
  : Returns a woods plot for comparisons of the timepoints samples
- [`woods_CI_plot_frac()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot_frac.md)
  : Returns a woods plot for comparisons of the timepoints samples
- [`deuteration_woods_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/deuteration_woods_timecourse.md)
  : Return woods plots for the timecourse
- [`deuteration_woods_timecourse_frac()`](https://mkajano.github.io/HDXBoxeR/reference/deuteration_woods_timecourse_frac.md)
  : Return woods plots for the timecourse
- [`average_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/average_timecourse.md)
  : Calculates average for time course data.
- [`sd_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timecourse.md)
  : Returns standard deviation for uptake data for timecourses.
- [`sd_timecourse_proc()`](https://mkajano.github.io/HDXBoxeR/reference/sd_timecourse_proc.md)
  : Returns standard deviation for percent deuteration data for
  timecourses.

## Legends and colors

- [`legend_heat_map()`](https://mkajano.github.io/HDXBoxeR/reference/legend_heat_map.md)
  : Legend for the heatmaps prep function.
- [`legend_heat_map_tc()`](https://mkajano.github.io/HDXBoxeR/reference/legend_heat_map_tc.md)
  : Legend for the heatmaps for timecourses.
- [`legend_heat_map_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/legend_heat_map_timecourse.md)
  : Legend for the heatmaps prep function for timecourses.
- [`legend_heat_map_tp()`](https://mkajano.github.io/HDXBoxeR/reference/legend_heat_map_tp.md)
  : Legend for the heatmaps.Extracts names from data.frame
- [`legend_heat_map_tp_proc()`](https://mkajano.github.io/HDXBoxeR/reference/legend_heat_map_tp_proc.md)
  : Legend for the heatmaps percent.Extracts names from data.frame
- [`legend_nm_bottom()`](https://mkajano.github.io/HDXBoxeR/reference/legend_nm_bottom.md)
  : Legend, bottom of the plots
- [`legend_raw_ave()`](https://mkajano.github.io/HDXBoxeR/reference/legend_raw_ave.md)
  : Legend for average plot.
- [`legend_raw_ave_proc()`](https://mkajano.github.io/HDXBoxeR/reference/legend_raw_ave_proc.md)
  : Preparatory function to draw legends for average procent
- [`legend_raw_ave_tc()`](https://mkajano.github.io/HDXBoxeR/reference/legend_raw_ave_tc.md)
  : Legend for average deuteration plot for timecourse.
- [`legend_sig_peptides()`](https://mkajano.github.io/HDXBoxeR/reference/legend_sig_peptides.md)
  : Legend for the significant peptides
- [`legend_states_PerD_bottom()`](https://mkajano.github.io/HDXBoxeR/reference/legend_states_PerD_bottom.md)
  : Legend, bottom of the plots
- [`legend_tc_bottom()`](https://mkajano.github.io/HDXBoxeR/reference/legend_tc_bottom.md)
  : Preparatory function returns legends for the timecourses.
- [`color_ranges_Blue_Red_heat_map()`](https://mkajano.github.io/HDXBoxeR/reference/color_ranges_Blue_Red_heat_map.md)
  : Returns color pallete from red to blue with number of colors for
  defined ranges
- [`color_ranges_Spectral()`](https://mkajano.github.io/HDXBoxeR/reference/color_ranges_Spectral.md)
  : Returns Spectral pallette with colors matching defined ranges
- [`reset_par()`](https://mkajano.github.io/HDXBoxeR/reference/reset_par.md)
  : Reset plotting window parameters to default

## Interactive plots and structures

- [`interactive_fissure_explorer()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_explorer.md)
  : Combined interactive explorer: volcano + fissure_vent + Woods +
  table

- [`interactive_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_vent.md)
  : Interactive fissure_vent plot, linked to the volcano/Woods/table
  views

- [`interactive_volcano()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_volcano.md)
  : Interactive volcano plot, linked to the fissure_vent/Woods/table
  views

- [`interactive_woods()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_woods.md)
  : Interactive Woods plot, linked to the volcano/fissure_vent/table
  views

- [`significant_peptides_table()`](https://mkajano.github.io/HDXBoxeR/reference/significant_peptides_table.md)
  : Interactive table of peptides, linked to the
  volcano/fissure_vent/Woods plots

- [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)
  : Map residue-level HDX data onto a 3D structure

- [`plot_hdx_structure_legend()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure_legend.md)
  :

  Legend for
  [`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)'s
  color scale

## PyMOL scripts

- [`pymol_script_average_residue()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_average_residue.md)
  : Writes a text files with pymol scripts to list significant residues.
- [`pymol_script_significant_peptide()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_peptide.md)
  : Writes a text files with pymol scripts to list significant peptides
- [`pymol_script_significant_peptide_proc()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_peptide_proc.md)
  : Writes a text files with pymol scripts to list significant peptides
- [`pymol_script_significant_residue()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_residue.md)
  : Writes a text files with pymol scripts to list significant residues.
- [`pymol_script_significant_residue_proc()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_residue_proc.md)
  : Writes a text files with pymol scripts to list significant residues.

## Reporting and export

- [`render_report()`](https://mkajano.github.io/HDXBoxeR/reference/render_report.md)
  : One-command HTML report of an HDX-MS comparison
- [`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md)
  : Summary table for HDX-MS reporting (Masson et al. 2019)
- [`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md)
  : Tidy table of peptide-level results
- [`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md)
  : Prepares output with HDX-MS data for publications
- [`verbose_timepoint_output()`](https://mkajano.github.io/HDXBoxeR/reference/verbose_timepoint_output.md)
  : Returns csv with averages from analysis for uptake file, standard
  deviation, p-values.
- [`verbose_timecourse_output()`](https://mkajano.github.io/HDXBoxeR/reference/verbose_timecourse_output.md)
  : Returns csv with averages from analysis for procent deuteration
  file, standard deviation for time courses.

## External tools

- [`extreme_input_gap()`](https://mkajano.github.io/HDXBoxeR/reference/extreme_input_gap.md)
  : Makes input for Extreme for bimodal analysis.
- [`extreme_input_undeut()`](https://mkajano.github.io/HDXBoxeR/reference/extreme_input_undeut.md)
  : Makes input for Extreme for bimodal analysis.
