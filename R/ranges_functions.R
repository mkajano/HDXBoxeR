#' Gives ranges for the averages
#'
#' Function used as internal function to get ranges in the function.
#'
#' @param df_ave average per residues
#' @param values_df data frame with values.
#' @return ranges per set
#' @export
ranges_function<-function(df_ave, values_df){
  values_df<-na.omit(values_df)
  message(paste("values range from ", round(range(values_df)[1],2), "to",
              round(range(values_df)[2],2), "%" ))

  lbs<-str_sub(colnames(df_ave[7:dim(df_ave)[2]]), start=4, end=-9)
  for ( i in 1:dim(values_df)[2]){
    lbs1=paste(lbs[1], lbs[i+1])
    message(paste(lbs1,"range", round(range(values_df[,i])[1],2),
                round(range(values_df[,i])[2],2)))
  }
}


#' Gives ranges for the averages for time course analysis
#'
#' Function used as internal function to get ranges in the function.
#'
#' @param df_ave average per residues
#' @param values_df data frame with values.
#' @return ranges per set
#' @export
ranges_function_tc<-function(df_ave, values_df){
  values_df<-na.omit(values_df)
  message(paste("values range from ", round(range(values_df)[1],2), "to",
              round(range(values_df)[2],2), "%" ))

  lbs<-str_sub(colnames(df_ave[7:dim(df_ave)[2]]), start=4, end=-9)
  for ( i in 1:dim(values_df)[2]){

    message(paste(lbs[i],"range", round(range(values_df[,i])[1],2),
                round(range(values_df[,i])[2],2)))
  }
}
