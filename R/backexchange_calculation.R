
######
#' Summary of backexchange summary
#'
#' Returns average and ranges of backexchange. Function calculates as: 1- (m100%-m0%)/N/Dfact.
#' m0% is the non-deuterated peptide centroid mass, m100% is the maximally labeled peptide centroid mass, N is
#' the theoretical number of backbone amides in the peptide and Dfrac is the fraction of D/H in the labeling buffer used.
#' Function requires undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD respectively.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export.
#' @param Dfact is the fraction of D/H in the labeling buffer used. Default set up to 0.85
#' @return Returns summary table for backexchange.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- backHX_calculations(filepath=file_nm, Dfact=0.85)
#' @export
backHX_calculations<- function(filepath, Dfact=0.85){

  a <- read.csv(file = filepath, header = F, skip = 1)
  nm <- read.csv(file = filepath, header = T, row.names = NULL,
                 nrows = 1)
  nm1 <- colnames(nm)

  if(length(unique(nm1=="row.names"))==2){
    nm1<-nm1[-which(nm1=="row.names")]
    nm1<-c(nm1, "row.names")}

  colnames(a) <- c(nm1)

  dif_col<- setdiff(c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                      "Charge", "Search.RT", "X..Deut", "Deut..", "Theor.Cent"), nm1)

  if (length(dif_col)!=0)
  {stop(paste("input file missing column(s): ", dif_col, sep=""))}


  a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
            "Charge", "Search.RT", "X..Deut", "Deut..", "Theor.Cent")]

  undeut<-a[which(a$Deut.Time == '0s'),]
  FD<-a[which(a$Deut.Time == 'FD'),]
#undeut<-na.omit(undeut)
#FD<-na.omit(FD)

tst<-merge(undeut, FD, by = c('Protein.State', 'Start','End', 'Sequence',
                   'Charge'))
bck_levels<-c()
for (state in unique(a$Protein.State)){
  temp1<-tst[which(tst$Protein.State ==state),]

dif_cent<-temp1$Theor.Cent.y-temp1$Theor.Cent.x
N<-nb_exch_deut(temp1)
backHX<-c((1-dif_cent/N/Dfact)*100)
bck_levels<-c(bck_levels,c(round(c(mean(backHX), range(backHX)), digits=2)))
#bck_levels<-c(bck_levels, paste(state, "state. Ave BackHX: ", round(mean(backHX), digits = 2), "% ; BackHX range: ",
#                  round(range(backHX)[1], digits = 2),"-", round(range(backHX)[2], digits = 2), "%"))
}
sum1<-data.frame(matrix(bck_levels,nrow=3))
names(sum1)<- unique(a$Protein.State)
row.names(sum1)<-c("<BackHX>", "BackHX Range From", "To")
return(t(sum1))}

