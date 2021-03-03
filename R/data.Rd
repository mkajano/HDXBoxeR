#' All_results_table.csv example input from HDXExaminer
#'
#' Info on how to prepare input from HDXExaminer
#' open your HDexaminer file
#' Under Peptides tab > Pool > right click peptide - delete peptides with no usable data
#' Tools > Options > Display > Deuteration Table - Display all columns (we did not select show low-confidence results)
#' Tools > Export > All results table
#'
#' @format A csv file:
#' \describe{
#'   \item{Protein State}{Name}
#'   \item{Deut Time}{Time in seconds}
#'   \item{Start}{number of first residue}
#'   \item{End}{last peptide residue number}
#'   \item{Sequence}{Amino Acid sequence}
#'   \item{Adduct}{}
#'   \item{Charge}{peptide charge}
#'   \item{Search RT}{Search Retention Time}
#'   \item{Actual RT}{Actual Retention Time}
#'   \item{# Spectra}{}
#'   \item{Search IMS}{}
#'   \item{IMS Range}{}
#'   \item{Max Inty}{}
#'   \item{Exp Cent}{}
#'   \item{Theor Cent}{}
#'   \item{Score}{}
#'   \item{Cent Diff}{}
#'   \item{# Deut}{Deuterons uptake}
#'   \item{% Deut}{Procent Deut}
#'   \item{Mark}{}
#'   \item{Confidence}{}
#' }
#' @source \url{Cook M, Delbecq SP, Schweppe TP, Guttman M, Klevit RE, Brzovic PS. The ubiquitin ligase SspH1 from Salmonella uses a modular and dynamic E3 domain to catalyze substrate ubiquitylation. J Biol Chem. 2019;294(3):783-793. doi:10.1074/jbc.RA118.004247}
"All_results_table"
