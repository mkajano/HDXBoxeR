## controls (fully deuterated "FD", in-exchange "0.00s") as one row per
## peptide and state, one column per replicate - built with the same engine
## as output_tc(); a single replicate is enough for a control
.hdx_control_table <- function(filepath, label, percent) {
  raw <- .hdx_raw(filepath)
  labels <- label
  if (label != "FD") {
    ## the in-exchange control can be written "0s", "0.00s", "0" ...
    tnum <- suppressWarnings(as.numeric(sub("s$", "", raw$Deut.Time)))
    zero <- unique(raw$Deut.Time[!is.na(tnum) & tnum == 0])
    if (length(zero)) labels <- zero
  }
  a <- raw[raw$Deut.Time %in% labels, .hdx_required, drop = FALSE]
  a$X..Deut <- suppressWarnings(as.numeric(a$X..Deut))
  a$Deut.. <- suppressWarnings(as.numeric(a$Deut..))
  a <- a[!is.na(if (percent) a$Deut.. else a$X..Deut), , drop = FALSE]
  if (nrow(a) == 0) stop("no '", paste(labels, collapse = "'/'"), "' rows with values in the input",
                         call. = FALSE)
  ## if several labels hold values, keep the one with the most rows
  label <- names(which.max(table(a$Deut.Time)))
  a <- a[a$Deut.Time == label, , drop = FALSE]
  out <- .hdx_wide(a, states = unique(a$Protein.State), times = label, min_replicates = 1,
                   percent = percent, layout = "tc")
  out[, c("Protein.State", "Start", "End", "Sequence", "Search.RT", "Charge",
          colnames(out)[7:ncol(out)])]
}

####################
#' Prepares output for HDX-MS Full deuteration data
#'
#' Returns the fully deuterated (`"FD"`) control: one row per peptide and
#' protein state, one column per FD replicate (uptake in Da).
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<-output_FD(file_nm)
#'
#' @export
output_FD<- function(filepath){
  .hdx_control_table(filepath, "FD", percent = FALSE)
}

####################
#' Prepares output for HDX-MS  Undeuterated sample data.
#'
#' Returns the in-exchange control (`Deut.Time` 0, e.g. `"0s"` or `"0.00s"`):
#' one row per peptide and protein state, one column per replicate (uptake in Da).
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return data frame with reorganized data where in columns is uptake data for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_UD(file_nm)
#' @export
output_UD<- function(filepath){
  .hdx_control_table(filepath, "0.00s", percent = FALSE)
}

####################
#' Prepares output for HDX-MS Full deuteration data  for procent deuteration.
#'
#' Returns the fully deuterated (`"FD"`) control as percent deuteration: one
#' row per peptide and protein state, one column per FD replicate.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return data frame with reorganized data where in columns is procent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_FD_proc(file_nm)
#' @export
output_FD_proc<- function(filepath){
  .hdx_control_table(filepath, "FD", percent = TRUE)
}

####################
#' Prepares output for HDX-MS Undeuterated data for procent deuteration.
#'
#' Returns the in-exchange control (`Deut.Time` 0, e.g. `"0s"`) as percent
#' deuteration: one row per peptide and protein state, one column per
#' replicate.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return data frame with reorganized data where in columns is procent deuteration for Protein States.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_UD_proc(file_nm)
#' @export
output_UD_proc<- function(filepath){
  .hdx_control_table(filepath, "0.00s", percent = TRUE)
}
