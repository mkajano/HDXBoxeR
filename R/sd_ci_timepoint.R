#' Calculate what would be the standard deviation based on the CI
#'
#' Calculate the CI to standard devaation, that can be used for the ellipsoid function
#' Calculate standard deviation from the CI using the equation as shown in:
#' SD=sqrt(replicates)*(uper_limit-lower_limit)/(2*t), t=qt(alpha/2, df) for
#' the two-sided (1-alpha) CI returned by t.test()
# https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals5.html
# https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
#'
#' @param df input from tp function
#' @param replicates number of replicates in sample. Default set to 3.
#' @param alpha what alpha (CI cutoff should be used)
#' @return data frame with standard deviation calculated from the CI
#' @export
sd_ci_timepoint<-function(df, replicates=3, alpha=0.01) {
  blocks <- .hdx_blocks(df, replicates)
  m1 <- .hdx_block_matrix(df, blocks[[1]])
  sd_one <- function(cols) {
    m2 <- .hdx_block_matrix(df, cols)
    vapply(seq_len(nrow(df)), function(i) {
      x1 <- m1[i, ][!is.na(m1[i, ])]; x2 <- m2[i, ][!is.na(m2[i, ])]
      if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
      n <- min(length(x1), length(x2))
      tvalue <- abs(qt(alpha/2, length(x1) + length(x2) - 2))
      tt <- tryCatch(t.test(x1, x2, conf.level = 1 - alpha), error = function(e) NULL)
      if (is.null(tt) || anyNA(tt$conf.int)) return(NA_real_)
      (tt$conf.int[2] - tt$conf.int[1]) / (2 * tvalue) * sqrt(n)
    }, numeric(1))
  }
  sd2 <- data.frame(lapply(blocks, sd_one))
  colnames(sd2) <- paste0("sd_ci_", names(blocks))
  sd2 <- data.frame(df[, 1:6], sd2)
  return(sd2)}
