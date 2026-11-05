# One-command HTML report of an HDX-MS comparison

Renders an HTML report with the analysis parameters, the Masson et al.
(2019) summary table
([`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md)),
data-quality checks
([`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md)
issues and the coverage map), significant-peptide counts and a
fissure-vent plot for every time point and state compared with the
reference, and - when given - region ANOVAs
([`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md))
and whole-time-course results
([`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md)).
The tidy results table
([`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md))
is written next to the report as csv, together with its parameters file.

## Usage

``` r
render_report(
  filepath,
  output_file = "hdx_report.html",
  states = NULL,
  times = NULL,
  reference = NULL,
  regions = NULL,
  alpha = 0.01,
  pvalue = 0.01,
  p_adjust = "none",
  percent = FALSE,
  metadata = list(),
  title = "HDX-MS report",
  Dfact = 0.85,
  protein_range = NULL
)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- output_file:

  path of the HTML report. Default `"hdx_report.html"` in the working
  directory.

- states, times, reference:

  states, time points and reference state (default: all states, all
  common time points, the first state).

- regions:

  optional named list of `c(start, end)` residue ranges for
  [`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md).

- alpha, pvalue, p_adjust, percent:

  analysis settings, as in
  [`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md).

- metadata:

  named list of extra summary-table items, as in
  [`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md).

- title:

  report title.

- Dfact:

  D2O fraction of the labeling buffer, for back-exchange. Default
  `0.85`.

- protein_range:

  residue range for the coverage, as in
  [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md).

## Value

the path of the report, invisibly.

## Details

Needs the rmarkdown and knitr packages (and pandoc, which comes with
RStudio).

## Examples

``` r
if (FALSE) { # \dontrun{
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
render_report(file_nm, output_file = file.path(tempdir(), "report.html"),
              regions = list(nterm = c(1, 100), cterm = c(101, 225)))
} # }
```
