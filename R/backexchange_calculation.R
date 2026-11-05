
######
#' Summary of backexchange summary
#'
#' Returns average and ranges of backexchange. Function calculates as: 1- (m100%-m0%)/N/Dfact.
#' m0% is the non-deuterated peptide centroid mass, m100% is the maximally labeled peptide centroid mass, N is
#' the theoretical number of backbone amides in the peptide and Dfrac is the fraction of D/H in the labeling buffer used.
#' Function requires undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD respectively.
#'
#' Masses are taken from the measured centroids: HDX Examiner's `Exp Cent`
#' column is an m/z value, so each centroid is multiplied by the peptide's
#' `Charge` to convert it to a mass, i.e.
#' `m100% - m0% = (Exp Cent[FD] - Exp Cent[0s]) * Charge` (in Da).
#' Earlier versions used the theoretical centroid (`Theor Cent`) without
#' the charge factor, which overstated back-exchange roughly in proportion to
#' the charge state. The returned table carries this as its `"method"`
#' attribute.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param Dfact is the fraction of D/H in the labeling buffer used. Default set up to 0.85
#' @return Returns summary table for backexchange.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- backHX_calculations(filepath=file_nm, Dfact=0.85)
#' @export
backHX_calculations<- function(filepath, Dfact=0.85){

  a <- .hdx_raw(filepath)
  a$Exp.Cent <- suppressWarnings(as.numeric(a$Exp.Cent))
  if (all(is.na(a$Exp.Cent))) {
    stop("back-exchange needs measured centroids (Exp Cent), which this input does not have",
         call. = FALSE)
  }


  a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
            "Charge", "Search.RT", "X..Deut", "Deut..", "Exp.Cent")]

  ## undeuterated rows are "0s" (HDXExaminer and read_dynamx()); "0.00s"
  ## is HDXExaminer's in-exchange control and is not used here
  a$Exp.Cent <- suppressWarnings(as.numeric(a$Exp.Cent))
  key <- c('Protein.State', 'Start', 'End', 'Sequence', 'Charge')
  ## one mean centroid per peptide, charge and state (replicates averaged),
  ## so the 0s/FD merge is one-to-one
  cent_mean <- function(d) {
    d <- d[!is.na(d$Exp.Cent), , drop = FALSE]
    if (nrow(d) == 0) return(data.frame(d[, key], Exp.Cent = numeric(0)))
    stats::aggregate(Exp.Cent ~ Protein.State + Start + End + Sequence + Charge, data = d, FUN = mean)
  }
  undeut <- cent_mean(a[a$Deut.Time == '0s', , drop = FALSE])
  FD <- cent_mean(a[a$Deut.Time == 'FD', , drop = FALSE])

tst<-merge(undeut, FD, by = key)
bck_levels<-c()
for (state in unique(a$Protein.State)){
  temp1<-tst[which(tst$Protein.State ==state),]

## Exp.Cent is m/z: multiply by charge to get the mass difference in Da
dif_cent<-(as.numeric(temp1$Exp.Cent.y)-as.numeric(temp1$Exp.Cent.x))*as.numeric(temp1$Charge)
N<-nb_exch_deut(temp1)
backHX<-c((1-dif_cent/N/Dfact)*100)
backHX <- backHX[is.finite(backHX)]
bck_levels<-c(bck_levels, if (length(backHX)) round(c(mean(backHX), range(backHX)), digits=2) else rep(NA_real_, 3))
}
sum1<-data.frame(matrix(bck_levels,nrow=3))
names(sum1)<- unique(a$Protein.State)
row.names(sum1)<-c("<BackHX>", "BackHX Range From", "To")
out<-t(sum1)
attr(out, "method")<-"measured centroid (Exp Cent) x Charge: 1 - (m_FD - m_0s)/N/Dfact"
return(out)}

