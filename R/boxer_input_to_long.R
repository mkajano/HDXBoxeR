
#' Prepares the long input of the input data
#'
#' The output is being used later for the analysis to compare different protein states.
#' input parameters match the output_tp parameters
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param replicates number of replicates to be used in analysis. The function takes number of replicates up to specified number. If no argument provided number maximal common number of replicates it used.
#' @param states function allows to choose what states should be used for analysis. Default all states are used.
#' @param times lists the deuteration times to be used in analysis. Default all states used.
#' @param seq_match Flag allows to choose if the peptide sequences should be matched between states. seq_match=FALSE signifies no sequence matching, seq_match=TRUE states that the sequences are matched between the sets.
#' @param csv Flag allowing saving the output as csv. With default csv="NA", data is not saved. If csv output is desided, provide output name.
#' @param percent Flag allowing to choose output as deteuration uptake (FALSE) or percent deuteration (TRUE). Default deuteration uptake.
#' @return data frame with reorganized in a long format
#' @importFrom tidyr gather
#' @export

boxer_input_to_long<-function(filepath, replicates, states, times, seq_match=FALSE, csv="NA", percent=FALSE){

  a5b<-output_tp(filepath, replicates, states, times, seq_match, csv, percent)

  v5 <- gather(a5b, Experiment, X..Deut,
               names(a5b)[7]:names(a5b)[length(names(a5b))], factor_key=TRUE)

  ## column names being split here are of the form "<state>_X..Deut_<n>"
  ## (percent=FALSE) or "<state>_Deut.._<n>" (percent=TRUE); pull the state
  ## and replicate number out with an explicit suffix match instead of
  ## splitting on a literal "X", which only exists in the uptake suffix.
  suffix_pattern <- "_(X\\.\\.Deut|Deut\\.\\.)_([0-9]+)$"
  v5$Protein.State <- sub(suffix_pattern, "", as.character(v5$Experiment))
  v5$Experiment <- paste("X..Deut", sub(paste0(".*", suffix_pattern), "\\2", as.character(v5$Experiment)), sep="_")

  v5$nb_prot<-nb_exch_deut(v5)
  if (percent==FALSE){
    v5$frac_deut<-as.numeric(v5$X..Deut)/v5$nb_prot
  } else {
    v5$frac_deut<-as.numeric(v5$X..Deut)/(100*v5$nb_prot)
  }

  return(v5)
}


#' Raw all_data input
#'
#' The output is being used later for the analysis to compare different protein states.
#' no initial matching of the files is made, raw data imported
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @return data frame with reorganized in a long format
#' @export
all_data_long_input<-function(filepath){
  b1<-.hdx_raw(filepath)
  bt<-b1[which(b1$Deut.Time != "0s" & b1$Deut.Time != "FD" ),]
  bt$nb_prot<-nb_exch_deut(bt)
  bt$frac_deut<-as.numeric(bt$X..Deut)/bt$nb_prot
  return(bt)
}


