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
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  sd_nm<-paste("sd_ci_", nm_root, sep="")

  sd1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {x1<-df[i,7:(7+replicates-1)];
    x2<-df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]

    tvalue=abs(qt(alpha/2, replicates*2-2))
    tt<-t.test(x1, x2, conf.level = c(1-alpha))
    sd.ci<-(tt$conf.int[2]-tt$conf.int[1])/(2*tvalue)*sqrt(replicates)
    sd1<-c(sd1, sd.ci)}
  }

  sd2<-data.frame(matrix(sd1, ncol=nb_sets , byrow = FALSE))
  colnames(sd2)<-sd_nm
  sd2<-data.frame(df[,1:6], sd2)
  return(sd2)}
