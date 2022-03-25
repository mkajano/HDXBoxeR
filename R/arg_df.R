#' Returns initially processed data.frame from the export from the HDXExaminer
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @return Data.frame for further processing
#' @export
arg_df<-function(filepath){
  a <- read.csv(file = filepath, header = F, skip = 1)
  nm <- read.csv(file = filepath, header = T, row.names = NULL, 
                 nrows = 1)
  nm1 <- colnames(nm)
  
  if(length(unique(nm1=="row.names"))==2){
    nm1<-nm1[-which(nm1=="row.names")]
    nm1<-c(nm1, "row.names")}
  
  colnames(a) <- c(nm1)
  
  dif_col<- setdiff(c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                      "Charge", "Search.RT", "X..Deut", "Deut.."), nm1)
  
  if (length(dif_col)!=0)
  {stop(paste("input file missing column(s): ", dif_col, sep=""))}
  
  
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
