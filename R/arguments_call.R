
#' Returns default arguments for the output_tp functions
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @return The default arguments to output_tp functions.
#' @export
arguments_call<-function(filepath){
  a <- read.csv(file = filepath, header = F, skip = 1)
  nm <- read.csv(file = filepath, header = T, row.names = NULL,
                 nrows = 1)
  nm1 <- colnames(nm)
  colnames(a) <- c(nm1)
  a <- a[, c(1:6)]
  if (all(a$Deut.Time == '0.00s')== FALSE & length(unique(a$Deut.Time == '0.00s'))==2){
    a<-a[-which(a$Deut.Time == c('0.00s')),]}
  if (all(a$Deut.Time == 'FD')== FALSE & length(unique(a$Deut.Time == 'FD'))==2){
    a<-a[-which(a$Deut.Time == c('FD')),]}
  if (all(a$Deut.Time == '0s')== FALSE & length(unique(a$Deut.Time == '0s'))==2){
    a<-a[-which(a$Deut.Time == c('0s')),]}
  a <- na.omit(a)

  states=unique(a$Protein.State)

  tm1<-unique(a[which(a$Protein.State==a$Protein.State[1]),which(colnames(a)=="Deut.Time")]); for (i in unique(a$Protein.State)){
    times=intersect(tm1, a[which(a$Protein.State==i),which(colnames(a)=="Deut.Time")])
  }


  len_rep<-c(); for (i in unique(a$Protein.State)){
    for (j in unique(a$Deut.Time)){

      len_rep<-c(len_rep, length(unique(a[which(a$Protein.State==i & a$Deut.Time==j),
                                          which(colnames(a)=="Experiment")])))
    }}; replicates=min(len_rep)

  return(list(states, times, replicates))}

