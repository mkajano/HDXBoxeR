#' Returns standard deviation for dataframe.
#'
#' Calculates the standard deviation of the replicates of every state (or
#' time, for [output_tc()] tables) for each peptide. Replicate blocks are
#' found from the column names, so states may have different numbers of
#' replicates and missing replicates (`NA`) are ignored; a peptide with
#' fewer than two values in a block gets `NA`.
#'
#' @param df output from functions output_tp or output_tc.
#' @param replicates only used for tables whose value columns are not named
#'   `<state>_..._<replicate>`: the number of replicates per block. Default 3.
#' @return Data.frame with standard deviation. The per-peptide number of
#'   replicates behind each SD is attached as attribute `"n"`.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' sd<-sd_timepoint(df=a, replicates=3)
#' @export
sd_timepoint<-function(df, replicates=3) {
  blocks <- .hdx_blocks(df, replicates)
  mats <- lapply(blocks, function(cols) .hdx_block_matrix(df, cols))
  sd2 <- data.frame(lapply(mats, .hdx_row_sd))
  colnames(sd2) <- paste0("sd_", names(blocks))
  sd2 <- data.frame(df[, 1:6], sd2)
  n <- data.frame(lapply(mats, .hdx_row_n)); colnames(n) <- names(blocks)
  rownames(n) <- rownames(sd2)
  attr(sd2, "n") <- n
  attr(sd2, "replicates") <- vapply(blocks, length, integer(1))
  class(sd2) <- c("hdx_sd", "data.frame")
  return(sd2)
}

#' Returns average value for either uptake of procent data.
#'
#' Calculates the average of the replicates of every state (or time, for
#' [output_tc()] tables) for each peptide, ignoring missing replicates.
#'
#' @param df output from functions output_tp or output_tc.
#' @param replicates only used for tables whose value columns are not named
#'   `<state>_..._<replicate>`: the number of replicates per block. Default 3.
#' @return Data.frame with average values
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' ave<-ave_timepoint(df=a)
#' @export
ave_timepoint<-function(df, replicates=3) {
  blocks <- .hdx_blocks(df, replicates)
  ave2 <- data.frame(lapply(blocks, function(cols) .hdx_row_mean(.hdx_block_matrix(df, cols))))
  colnames(ave2) <- paste0("av_", names(blocks))
  ave2 <- data.frame(df[, 1:6], ave2)
  attr(ave2, "replicates") <- vapply(blocks, length, integer(1))
  return(ave2)
}

#' Calculation of pvalue between first protein state and any other state from all_states file
#'
#' Welch t-test of every state (or time) against the first one in `df`, per
#' peptide. Missing replicates are ignored; the p-value is `NA` when either
#' side has fewer than two values or the test is undefined (e.g. all
#' replicates identical).
#'
#' @param df output from functions output_tp or output_tc.
#' @param replicates only used for tables whose value columns are not named
#'   `<state>_..._<replicate>`: the number of replicates per block. Default 3.
#' @return Data.frame with p-values
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' pv<-pv_timepoint(df=a)
#' @export
pv_timepoint<-function(df, replicates=3) {
  blocks <- .hdx_blocks(df, replicates)
  m1 <- .hdx_block_matrix(df, blocks[[1]])
  pv2 <- data.frame(lapply(blocks, function(cols) .hdx_row_welch(m1, .hdx_block_matrix(df, cols))))
  colnames(pv2) <- paste0("pv_", names(blocks))
  pv2 <- data.frame(df[, 1:6], pv2)
  return(pv2)}


#' Returns data frame with difference of averages between State1 and other states provided.
#'
#' Returns average difference data.frame.
#' Sets are compared to the first state in the input file.
#' If other order of the sets is required use
#' Default for the number of replicates is 3.
#'
#' @param df output from functions output_tp, output_tp_proc,  output_tp_states or output_tp_proc_states.
#' @return Data.frame with difference values btw control and other protein states.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' pv<-pv_timepoint(df=a) ##if number of replicates is equal 3
#' pv1<-pv_timepoint(df=a, replicates=3) ##if number of replicates is equal 4
#' #b<-output_tp_states(file_nm, states=c("4EHP", "State2", "State3" ))
#' #pv_states<-pv_timepoint(df=b) ### here means of State4, will be compared to State2 and State4
#' @export
dif_ave<-function(df){
  da1<-data.frame(df[,1:6],df[,7:dim(df)[2]]-df[,7])
  return(da1)
}
