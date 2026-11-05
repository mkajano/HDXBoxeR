#' Returns initially processed data.frame from the export from the HDXExaminer
#'
#' Function used as internal function
#'
#' @param filepath input file location, or a table from [read_hdx()]
#' @return Data.frame for further processing
#' @export
arg_df<-function(filepath){
  a <- .hdx_raw(filepath)
  a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
            "Charge", "Search.RT", "X..Deut", "Deut..")]
  if (all(a$Deut.Time == "0.00s") == FALSE & length(unique(a$Deut.Time ==
                                                           "0.00s")) == 2) {
    a <- a[-which(a$Deut.Time == c("0.00s")), ]  }
  if (all(a$Deut.Time == "FD") == FALSE & length(unique(a$Deut.Time ==
                                                        "FD")) == 2) {
    a <- a[-which(a$Deut.Time == c("FD")), ]}
  if (all(a$Deut.Time == "0s") == FALSE & length(unique(a$Deut.Time ==
                                                        "0s")) == 2) {
    a <- a[-which(a$Deut.Time == c("0s")), ] }
  a <- na.omit(a)


  return(a)
}
