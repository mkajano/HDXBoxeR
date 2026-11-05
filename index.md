# HDXBoxeR

![figure_bioinformatics2](https://github.com/mkajano/HDXBoxeR/assets/11892713/2b16a28f-dcfa-40d6-af04-08704dfd3c32)

figure_bioinformatics2

HDXBoxeR analyses hydrogen-deuterium exchange mass spectrometry (HDX-MS)
data: from the exported peptide table to statistically tested
differences between protein states, publication tables, figures and
structure-coloring scripts.

## What it does

- **Reads** HDXExaminer “All results” tables and Waters DynamX
  cluster/state exports
  ([`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)),
  or any long table
  ([`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)).
  States may have different numbers of replicates, and missing
  replicates are handled.
- **Tests** differences between states with one rule throughout - Welch
  t-test p-value AND the pooled critical interval of Hageman &
  Weis (2019) - with optional false-discovery-rate correction
  (`p_adjust = "BH"`) and length-normalized values
  (`percent = "per_residue"`).
- **Compares** two states peptide by peptide
  ([`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
  volcano, Woods, heat maps, robot plots), three or more states
  ([`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)),
  whole regions
  ([`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md),
  a region ANOVA with peptides as their own controls), whole time
  courses
  ([`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md)),
  states with different peptide maps
  ([`compare_states_unmatched()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_unmatched.md)),
  and finds contiguous affected regions
  ([`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)).
- **Checks data quality**: coverage, redundancy, back-exchange, and
  peptide-level flags
  ([`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md),
  [`plot_peptide_coverage_map()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_coverage_map.md)).
- **Normalizes** to the fully deuterated control (`percent = "FD"`).
- **Reports**: the summary table recommended by Masson et al. (2019)
  ([`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md)),
  a tidy results table with its parameters
  ([`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md)),
  publication-format data
  ([`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md)),
  and a one-command HTML report
  ([`render_report()`](https://mkajano.github.io/HDXBoxeR/reference/render_report.md)).
- **Visualizes on structures**: PyMOL scripts and interactive NGL views
  ([`plot_hdx_structure()`](https://mkajano.github.io/HDXBoxeR/reference/plot_hdx_structure.md)),
  plus linked interactive plots
  ([`interactive_fissure_explorer()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_explorer.md)).

## Installation

From CRAN:

``` r

install.packages("HDXBoxeR")
```

Development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("mkajano/HDXBoxeR")
```

## Quick start

``` r

library(HDXBoxeR)
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")

d <- read_hdx(file_nm)                 # HDXExaminer or DynamX; or as_hdx_long() for other software
hdx_summary_table(d)                   # Masson et al. 2019 summary table
peptide_qc(d)                          # peptide-level quality flags

a <- output_tp(d)                      # one row per peptide and time, one column per replicate
plot_fissure_vent(a, times = "60.00s", state1 = "Unbound", state2 = "bound")
res <- hdx_results_table(d, p_adjust = "BH")   # every peptide, time and pair, with the call

compare_regions(d, regions = list(nterm = c(1, 100), cterm = c(101, 225)))

render_report(d, output_file = "hdx_report.html",
              regions = list(nterm = c(1, 100), cterm = c(101, 225)))
```

Two example data sets come with the package: `All_results_table.csv`
(two states) and `three_states.csv` (wild type and two point mutants,
uneven replicates; Woods et al. 2025, PNAS) - see
[`?All_results_table`](https://mkajano.github.io/HDXBoxeR/reference/All_results_table.md)
and
[`?three_states`](https://mkajano.github.io/HDXBoxeR/reference/three_states.md).

## Vignettes

- [`vignette("HDXBoxeR")`](https://mkajano.github.io/HDXBoxeR/articles/HDXBoxeR.md) -
  the main workflow: data preparation, statistics, plots, PyMOL scripts,
  summary tables
- [`vignette("input_formats")`](https://mkajano.github.io/HDXBoxeR/articles/input_formats.md) -
  HDXExaminer, DynamX and other software; replicates; FD normalization
- [`vignette("fissure_vent")`](https://mkajano.github.io/HDXBoxeR/articles/fissure_vent.md) -
  the fissure-vent plot
- [`vignette("multistate")`](https://mkajano.github.io/HDXBoxeR/articles/multistate.md) -
  three or more states, whole regions and time courses
- [`vignette("data_quality")`](https://mkajano.github.io/HDXBoxeR/articles/data_quality.md) -
  coverage, redundancy and back-exchange
- [`vignette("interactive_explorer")`](https://mkajano.github.io/HDXBoxeR/articles/interactive_explorer.md) -
  linked interactive plots
- [`vignette("hdx_structure")`](https://mkajano.github.io/HDXBoxeR/articles/hdx_structure.md) -
  mapping results onto a structure

## Citation

Janowska MK, Reiter K, Magala P, Guttman M, Klevit RE (2024). HDXBoxeR:
an R package for statistical analysis and visualization of multiple
hydrogen-deuterium exchange mass-spectrometry datasets of different
protein states. *Bioinformatics* 40(8):btae479.
<doi:10.1093/bioinformatics/btae479>

See `citation("HDXBoxeR")`.

## License

GPL (\>= 2)
