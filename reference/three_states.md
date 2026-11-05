# Example input: three protein states with point mutations

A second example HDXExaminer "All results" export, for multi-state
comparisons: a wild-type protein and two point mutants, `WT_pH7.5`,
`D109H_pH7.5` and `R120G_ph7.5`. It is a file, not an R object - get its
path with
`system.file("extdata", "three_states.csv", package = "HDXBoxeR")`.

## Format

A csv file with the same columns as
[All_results_table](https://mkajano.github.io/HDXBoxeR/reference/All_results_table.md)
(without `Adduct`, `Search IMS` and `IMS Range`).

## Source

Woods et al. (2025) Proc Natl Acad Sci USA.
[doi:10.1073/pnas.2425061122](https://doi.org/10.1073/pnas.2425061122)

## Details

Labeling times 4 s, 60 s and 1800 s, plus undeuterated (`0s`) and fully
deuterated (`FD`) controls (a subset of the full data set). It shows
features of real data the first example lacks:

- uneven replicates: 3 or 4 per state and time, and only 2 undeuterated
  controls for R120G;

- peptides covering residue 109 or 120 differ in sequence between the
  states (D109H, R120G), so they are matched only with
  `seq_match = TRUE` / `by = "position"`.
