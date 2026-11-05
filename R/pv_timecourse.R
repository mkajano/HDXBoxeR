#' pvalue calculation between two sets of the data at certain timepoint
#'
#' Preparatory function for calculation of pvalue between sets.
#'
#' @param df_c dataframe of control
#' @param df_v dataframe for variant
#' @param replicates number of replicates. Default set to 3.
#' @param p_adjust multiple-testing correction passed to [stats::p.adjust()]:
#'   `"none"` (default) or e.g. `"BH"`, applied across the peptides of the
#'   comparison, separately for each time point.
#' @return pvalue comparisons between two sets.
#' @keywords internal
#' @export
pv_timecourse<-function(df_c,df_v, replicates=3, p_adjust="none") {
  keys <- c('Start','End', 'Sequence', 'Search.RT', 'Charge')
  blocks_c <- .hdx_blocks(df_c, replicates)
  blocks_v <- .hdx_blocks(df_v, replicates)
  ## match peptides between the two tables by row position, keeping merge()'s
  ## row order and its first six output columns as before
  ic <- data.frame(df_c[, keys], .row_c = seq_len(nrow(df_c)))
  iv <- data.frame(df_v[, keys], .row_v = seq_len(nrow(df_v)))
  combined_df <- merge(df_c, df_v, by = keys)
  idx <- merge(ic, iv, by = keys)
  pv2 <- data.frame(lapply(seq_along(blocks_c), function(j) {
    .hdx_row_welch(.hdx_block_matrix(df_c, blocks_c[[j]])[idx$.row_c, , drop = FALSE],
                   .hdx_block_matrix(df_v, blocks_v[[j]])[idx$.row_v, , drop = FALSE])
  }))
  colnames(pv2) <- paste0("pv_", names(blocks_c))
  pv2<-data.frame(combined_df[,1:6], pv2)
  return(.hdx_adjust_pv(pv2, p_adjust, skip_first = FALSE))
}
