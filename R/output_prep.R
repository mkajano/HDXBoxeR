#' Prepares output with HDX-MS data for publications
#'
#' Format prepared based of example from:
#' Masson, G.R., Burke, J.E., Ahn, N.G. et al. Recommendations for performing, interpreting and reporting hydrogen deuterium exchange mass spectrometry (HDX-MS) experiments. Nat Methods 16, 595–602 (2019). https://doi.org/10.1038/s41592-019-0459-y
#' It generates csv file in format ready for publication of the data.
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param output_name Name of output file. It has to be csv file
#' @param replicates optional maximum number of replicates per state and time. Default `NULL` uses every replicate; states may differ.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param percent return either uptake or percent deuteration, default=FALSE, return uptake
#' @return Writes the table to `output_name` and returns it invisibly.
#'   `Mean.Peptide_Mass_[Da]` is the mean measured neutral peptide mass,
#'   Charge x (centroid m/z - proton mass); earlier versions reported the
#'   centroid m/z under this name.
#' @importFrom grDevices col2rgb colorRampPalette rgb
#' @importFrom graphics abline axis box boxplot legend mtext par plot points polygon rect text arrows
#' @importFrom stats df na.omit qt sd t.test
#' @importFrom utils head read.csv tail type.convert write.csv write.table
#' @importFrom dplyr arrange %>%
#' @importFrom stringr str_sub str_count
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr separate
#' @importFrom methods rbind2
#' @importFrom wrapr orderv
#' @importFrom DescTools DrawEllipse
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' output_prep(filepath=file_nm, output_name=tempfile())
#' }
#' @export
output_prep<-function(filepath, output_name, states, replicates = NULL, times, percent=FALSE){
  if(missing(states)) { states=arguments_call1(filepath); .hdx_msg("Protein.States used:", states)}
  if(missing(times)) { times=arguments_call2(filepath, states); .hdx_msg("Deut.times used:", times)}

  ctl <- arg_UN_FD(filepath)
  undeut <- ctl[[1]]
  FD <- ctl[[2]]

  ## one row per measurement; mean and SD over the replicates (Experiments)
  ## of each state, time and peptide - states may have different numbers of
  ## replicates, and a missing centroid does not remove the uptake value
  a <- .hdx_raw(filepath)
  a <- a[a$Protein.State %in% states & a$Deut.Time %in% times, , drop = FALSE]
  a$val <- suppressWarnings(as.numeric(if (percent == FALSE) a$X..Deut else a$Deut..))
  a <- a[!is.na(a$val), , drop = FALSE]
  if (!is.null(replicates)) {
    keep <- unlist(lapply(split(seq_len(nrow(a)), list(a$Protein.State, a$Deut.Time), drop = TRUE), function(i) {
      i[a$Experiment[i] %in% utils::head(unique(a$Experiment[i]), replicates)]
    }))
    a <- a[sort(keep), , drop = FALSE]
  }
  n_rep <- tapply(a$Experiment, a$Protein.State, function(e) length(unique(e)))
  .hdx_msg("Number of replicates used:", paste0(names(n_rep), "=", n_rep))

  ## Exp.Cent is a centroid m/z: the peptide mass column reports the neutral
  ## mass, Charge * (m/z - proton mass), in Da
  to_mass <- function(mz, z) as.numeric(z) * (suppressWarnings(as.numeric(mz)) - .proton_mass)
  a$mass <- to_mass(a$Exp.Cent, a$Charge)
  g <- interaction(a$Protein.State, a$Deut.Time, a$Start, a$End, a$Sequence, a$Charge, drop = TRUE, sep = "\r")
  first <- !duplicated(g)
  m <- function(v) { x <- tapply(v, g, mean, na.rm = TRUE); x[is.nan(x)] <- NA; as.numeric(x[as.character(g[first])]) }
  s <- function(v) as.numeric(tapply(v, g, function(x) if (sum(!is.na(x)) < 2) NA_real_ else stats::sd(x, na.rm = TRUE))[as.character(g[first])])
  rt <- as.numeric(tapply(suppressWarnings(as.numeric(a$Search.RT)), g, stats::median, na.rm = TRUE)[as.character(g[first])])
  b <- a[first, c("Protein.State", "Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge")]
  b$Search.RT <- rt
  out <- data.frame(b, m(a$mass), s(a$mass), round(m(a$val), digits = 4), round(s(a$val), digits = 4))

  ## FD's "mean" column is not a real average: undeut is defined as 0
  ## exchange in both modes; FD uses its own measured X..Deut (Da) in uptake
  ## mode, or the conventional 100% in percent mode. SD stays 0.
  FD[,9] <- to_mass(FD[,9], FD$Charge)
  undeut[,9] <- to_mass(undeut[,9], undeut$Charge)
  FD2<-data.frame(FD[,c(1,2,4:6,8,7,9)], rep(0,times=dim(FD)[1]),
                  if (percent==FALSE) as.numeric(FD[,10]) else rep(100,times=dim(FD)[1]),
                  rep(0,times=dim(FD)[1]))
  undeut2<-data.frame(undeut[,c(1,2,4:6,8,7,9)], rep(0,times=dim(undeut)[1]),  rep(0,times=dim(undeut)[1]), rep(0,times=dim(undeut)[1]))
  FD2 <- FD2[FD2[[1]] %in% states, , drop = FALSE]
  undeut2 <- undeut2[undeut2[[1]] %in% states, , drop = FALSE]

  nm_final<-c("Protein.State", "Deut.Time", "Start", "End", "Sequence", "Retention_Time_[min]","Charge",
              "Mean.Peptide_Mass_[Da]", "st.dev_Peptide_Mass",
              if (percent == FALSE) c("Mean.Deut.Uptake_[Da]", "st.dev_Deut_Uptake")
              else c("Mean.Percent_Deut_[%]", "st.dev_Percent_Deut"))
  colnames(FD2)<-nm_final
  colnames(undeut2)<-nm_final
  colnames(out)<-nm_final

  out<-rbind(out, FD2)
  out<-rbind(out, undeut2)
  out<-out[,c(1:8,10:11)]


  out<-(arrange(out, Start, End, Protein.State, Charge))
  ##find a way to order correctly assign negative time to zero and
  ##very long time for FD and order based of these columns
  ##if bp does not work use out
  ord=out$Deut.Time
  ord=type.convert(ord, as.is = TRUE)
  ord[ord=="0s"]<- c("-5.0")
  ord[ord== "FD"]<- c("10000000000000.0s")
  ord=as.numeric(str_sub(ord, end=-2))

  bp<-data.frame(out, ord, check.names = FALSE)
  bp<-arrange(bp, Start, End, Protein.State, Charge, ord)
  bp<-bp[,-dim(bp)[2]]
  rownames(bp) <- NULL
  ###write output
  write.csv(bp, output_name, row.names = FALSE)
  invisible(bp)}




#' Returns initially processed data.frame from the export from the HDXExaminer
#'
#' Function used as internal function
#'
#' @param filepath input file location
#' @return Data.frame for further processing
#' @keywords internal
#' @export
arg_UN_FD<-function(filepath){
  a <- .hdx_raw(filepath)

  undeut <- a[a$Deut.Time=="0s",c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                                  "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]
  FD <- a[a$Deut.Time=="FD",c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                              "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]


  a <- a[,c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
            "Charge", "Search.RT","Exp.Cent", "X..Deut", "Deut..")]
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


  return(list(undeut, FD, a))
}


