#' Peptides significantly different between states
#'
#' Returns list of the dataframes with significantly different peptides
#' @param df differences in averages data.frame calculated using diff_ave function
#' @param alpha critical interval, here is multiple sets are using maximum CI is used.
#' @param replicates number of replicates
#' @param pvalue p-value cutoff here set up to 0.01
#' @param pv_cutoff deprecated, use `pvalue`.
#' @param p_adjust multiple-testing correction of the Welch p-values, passed to [stats::p.adjust()]: `"none"` (default, raw p-values) or e.g. `"BH"`. A peptide is significant when its (adjusted) p-value is below the cutoff AND its difference exceeds the critical interval.
#' @return list of the dataframes
#' @export
diff_peptides_prep<-function(df,replicates=3, pvalue=0.01, alpha=0.01, p_adjust="none", pv_cutoff){
  if (!missing(pv_cutoff)) { .hdx_deprecated("pv_cutoff", "pvalue"); pvalue <- pv_cutoff; if (missing(alpha)) alpha <- pv_cutoff }

  av1<-ave_timepoint(df, replicates)
  df1<-dif_ave(av1)
  pv<-pv_timepoint(df, replicates, p_adjust = p_adjust)
  s1<-sd_timepoint(df,  replicates)

  nb1=0
  peptide_list<-list()
  for ( k in(unique(df1$Deut.Time))){


    nb1=nb1+1
    nb2=0
    peptide_list_state<-list()
    for ( i in 8:dim(df1)[2]){

      nb2=nb2+1
      CI<-CI_tp(s1[s1$Deut.Time==k,c(1:7,i)], replicates, alpha)
      nm_of_state<-c(paste(.hdx_strip_suffix(colnames(df1)[7]),
                     " - ",
                     .hdx_strip_suffix(colnames(df1)[i])))

      abs.a<-abs(df1[df1$Deut.Time==k,i])>CI
      cl1<-pv[pv$Deut.Time==k,i]< pvalue
      df_sig_pep<-c(abs.a*cl1)

      ind_sig<-which(df_sig_pep==1)

      p_values_peptides<- pv[pv$Deut.Time==k,i][ind_sig]
      delta_uptake_from_control<-df1[df1$Deut.Time==k,i][ind_sig]

      df_sig_peptides<-data.frame(df1[df1$Deut.Time==k,1:6][ind_sig,],
                                  p_values_peptides,
                                  delta_uptake_from_control)




      peptide_list_state[[nb2]]<- list(nm_of_state, k, CI, df_sig_peptides)
      names(peptide_list_state[[nb2]]) <- c("names_of_states",
                                   "Deut_time", "Critical Interval",
                                   "significantly diff peptides")

    }
    peptide_list[[nb1]]<-peptide_list_state
  }
    return(peptide_list)
  }



