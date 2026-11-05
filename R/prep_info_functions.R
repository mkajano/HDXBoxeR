######
#' Provides summary table for all data.sets.
#'
#' Returns a data frame summarizing each protein state: time points, number
#' of replicates, number of peptides, sequence coverage, average peptide
#' length and redundancy.
#'
#' Peptides are those found at every time point of the state with at least
#' two replicates (as in [output_tp()]); each peptide sequence counts once
#' for coverage, length and redundancy (charge states are not counted
#' twice).
#' \itemize{
#'   \item coverage: percent of residues in `protein_range` covered by at
#'     least one peptide;
#'   \item redundancy: average number of peptides covering a covered residue;
#'   \item `# replicates`: replicates per time point (a range when they
#'     differ).
#' }
#'
#' @param filepath filepath to the input file. Input file is All_results table from HDX_Examiner, where all the fields are marked for export. Can also be a DynamX cluster/state export, or a table from [read_hdx()] / [as_hdx_long()].
#' @param protein_range residue range for the coverage, `c(first, last)`.
#'   Default: from the first to the last residue covered by any peptide in
#'   the data (the same range for every state).
#' @return Returns summary table.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- general_info(file_nm)
#' general_info(file_nm, protein_range = c(1, 225))
#' @export
general_info<- function(filepath, protein_range = NULL){
  raw <- .hdx_raw(filepath)
  a <- arg_df(raw)
  states <- unique(a$Protein.State)
  if (is.null(protein_range)) protein_range <- c(min(a$Start), max(a$End))
  span <- protein_range[1]:protein_range[2]

  rows <- lapply(states, function(state) {
    t_state <- unique(a$Deut.Time[a$Protein.State == state])
    t_state <- t_state[order(suppressWarnings(as.numeric(sub("s$", "", t_state))))]
    w <- suppressMessages(output_tp(raw, states = state, times = t_state))
    pep <- unique(w[, c("Start", "End", "Sequence")])
    n_rep <- tapply(a$Experiment[a$Protein.State == state], a$Deut.Time[a$Protein.State == state],
                    function(e) length(unique(e)))
    cov <- integer(max(c(span, pep$End)))
    for (i in seq_len(nrow(pep))) cov[pep$Start[i]:pep$End[i]] <- cov[pep$Start[i]:pep$End[i]] + 1L
    cov_span <- cov[span]
    data.frame(
      Protein.State = state,
      Timepoints = paste(t_state, collapse = " "),
      `# replicates` = if (min(n_rep) == max(n_rep)) as.character(min(n_rep)) else paste0(min(n_rep), "-", max(n_rep)),
      `# peptides` = nrow(unique(w[, c("Start", "End", "Sequence", "Charge")])),
      `peptide coverage %` = round(100 * mean(cov_span > 0), 2),
      `<Peptide Length>` = round(mean(pep$End - pep$Start + 1), 2),
      `<Redundancy>` = round(sum(cov_span) / sum(cov_span > 0), 2),
      check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- states
  attr(out, "protein_range") <- protein_range
  out
}
