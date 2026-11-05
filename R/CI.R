#' Global confidence interval treshold from experimental standard deviation for 2 samples.
#'
#' Calculation of global confidence interval using approach by:
#' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#' calculations for alpha 0.99
#'
#' Unequal or missing replicates: each state's pooled SD weights every
#' peptide's SD by its degrees of freedom (n_i - 1), the standard error of
#' the difference is `sqrt(sp1^2/n1 + sp2^2/n2)` with n1, n2 the average
#' number of replicates per peptide in each state, and the t-value uses
#' df = n1 + n2 - 2. With equal, complete replicates this is exactly the
#' Hageman & Weis formula.
#'
#' @param s1 standard deviation from one sample
#' @param s2 standard deviation from seconda sample
#' @param replicates number of replicates. A single number, or two numbers
#'   (state 1, state 2) when the states have different replicate counts.
#'   Default set to 3.
#' @param alpha probability to calculate t-value in Welch test. Default 0.01
#' @param n1,n2 optional per-peptide numbers of replicates behind `s1` and
#'   `s2` (e.g. from attribute `"n"` of [sd_timepoint()]). Default: every
#'   peptide has `replicates` values.
#' @return treshold for determining significance.
#' @examples
#' sd1<-data.frame(c(0.1, 0.12, 0.13, 0.09, 0.11, 0.10))
#' sd2<-data.frame(c(0.18, 0.11, 0.13, 0.08, 0.11, 0.06))
#' CI_2pts(s1=sd1, s2=sd2, replicates=3, alpha=0.01)
#' @export
CI_2pts<-function(s1, s2, replicates=3, alpha=0.01, n1=NULL, n2=NULL){
  s1 <- as.numeric(unlist(s1)); s2 <- as.numeric(unlist(s2))
  if (is.null(n1)) n1 <- rep(replicates[1], length(s1))
  if (is.null(n2)) n2 <- rep(replicates[length(replicates)], length(s2))
  sp1 <- .hdx_pooled_sd(s1, n1)
  sp2 <- .hdx_pooled_sd(s2, n2)
  N1 <- mean(n1[!is.na(s1) & n1 >= 2]); N2 <- mean(n2[!is.na(s2) & n2 >= 2])
  tvalue <- abs(qt(alpha, N1 + N2 - 2))
  spa <- sqrt(sp1^2 / N1 + sp2^2 / N2)
  CI <- spa * tvalue
  return(CI)
}

#' Global confidence interval treshold from experimental standard deviation
#'
#' Calculation of global confidence interval using approach by for all protein states compared to first state in the data.frame.
#' Reliable Identification of Significant Differences in Differential Hydrogen Exchange-Mass Spectrometry Measurements Using a Hybrid Significance Testing Approach
#' Tyler S. Hageman and David D. Weis
#' Analytical Chemistry 2019 91 (13), 8008-8016 DOI: 10.1021/acs.analchem.9b01325
#'
#'
#' @param df standard deviation dataframe.
#' @param replicates number of replicates. Default set to 3.
#' @param alpha significance level. Set as default to 0.01
#' @return treshold for determining significance.
#' @examples
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tc(file_nm, seq_match=FALSE)
#' sd<-sd_timepoint(df=a, replicates=3)
#' CI_tp(df=sd, replicates=3, alpha=0.01 )
#' CI_tp(sd)
#' @export
CI_tp<-function(df, replicates=3, alpha=0.01 ){
  CI_all<-c()
  for ( i in 8:dim(df)[2]){
    CI_all<-c(CI_all, CI_2pts(df[,7], df[,i], replicates, alpha,
                              n1 = .hdx_n_for(df, colnames(df)[7]),
                              n2 = .hdx_n_for(df, colnames(df)[i])))}
  return(CI_all)
}
