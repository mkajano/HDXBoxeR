#' Example input: All_results_table.csv from HDXExaminer
#'
#' The package ships one example HDXExaminer "All results" export, used by
#' the examples, tests and vignettes. It is a file, not an R object - get its
#' path with
#' `system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")`.
#' It holds two protein states ("Unbound", "bound"), four labeling times
#' (3 s, 60 s, 1800 s, 72000 s) plus undeuterated (`0s`) and fully
#' deuterated (`FD`) controls, and three replicates per state and time.
#'
#' How to prepare this export from your own HDXExaminer project:
#' \enumerate{
#'   \item Peptides tab > Pool > right-click a peptide - delete peptides with no usable data.
#'   \item Tools > Options > Display > Deuteration Table - display all columns
#'     (low-confidence results were not shown for this file).
#'   \item Tools > Export > All results table.
#' }
#'
#' @format A csv file with one row per peptide, charge state, protein state,
#'   labeling time and replicate:
#' \describe{
#'   \item{Protein State}{protein state name}
#'   \item{Deut Time}{labeling time, e.g. `"3.00s"`; `"0s"` undeuterated, `"FD"` fully deuterated}
#'   \item{Experiment}{replicate (experiment) name}
#'   \item{Start}{first residue number of the peptide}
#'   \item{End}{last residue number of the peptide}
#'   \item{Sequence}{amino acid sequence}
#'   \item{Adduct}{adduct ion}
#'   \item{Charge}{peptide charge state}
#'   \item{Search RT}{search retention time}
#'   \item{Actual RT}{actual retention time range}
#'   \item{# Spectra}{number of spectra}
#'   \item{Search IMS}{search ion-mobility value}
#'   \item{IMS Range}{ion-mobility range}
#'   \item{Peak Width}{peak width}
#'   \item{m/z Shift}{m/z shift}
#'   \item{Max Inty}{maximum intensity}
#'   \item{Exp Cent}{measured centroid (m/z)}
#'   \item{Theor Cent}{theoretical centroid (m/z)}
#'   \item{Score}{fit score}
#'   \item{Cent Diff}{centroid difference}
#'   \item{# Deut}{deuterium uptake (Da)}
#'   \item{Deut %}{percent deuteration}
#'   \item{Mark}{user mark}
#'   \item{Confidence}{HDXExaminer confidence level}
#' }
#' @source Cook M, Delbecq SP, Schweppe TP, Guttman M, Klevit RE, Brzovic PS.
#'   The ubiquitin ligase SspH1 from Salmonella uses a modular and dynamic E3
#'   domain to catalyze substrate ubiquitylation. J Biol Chem.
#'   2019;294(3):783-793. \doi{10.1074/jbc.RA118.004247}
#' @name All_results_table
#' @keywords datasets
NULL

#' Example input: three protein states with point mutations
#'
#' A second example HDXExaminer "All results" export, for multi-state
#' comparisons: a wild-type protein and two point mutants, `WT_pH7.5`,
#' `D109H_pH7.5` and `R120G_ph7.5`. It is a file, not an R object - get its
#' path with
#' `system.file("extdata", "three_states.csv", package = "HDXBoxeR")`.
#'
#' Labeling times 4 s, 60 s and 1800 s, plus undeuterated (`0s`) and fully
#' deuterated (`FD`) controls (a subset of the full data set). It shows
#' features of real data the first example lacks:
#' \itemize{
#'   \item uneven replicates: 3 or 4 per state and time, and only 2
#'     undeuterated controls for R120G;
#'   \item peptides covering residue 109 or 120 differ in sequence between
#'     the states (D109H, R120G), so they are matched only with
#'     `seq_match = TRUE` / `by = "position"`.
#' }
#'
#' @format A csv file with the same columns as [All_results_table]
#'   (without `Adduct`, `Search IMS` and `IMS Range`).
#' @source Woods et al. (2025) Proc Natl Acad Sci USA.
#'   \doi{10.1073/pnas.2425061122}
#' @name three_states
#' @keywords datasets
NULL
