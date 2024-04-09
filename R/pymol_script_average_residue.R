#' Writes a text files with pymol scripts to list significant residues.
#'
#' Function write a script that can be used in pymol to color structure.
#' Number of colors and corresponding to them ranges can be defined by user.
#' Residues are being colored by average uptake values from the significant peptides per residues.
#'
#' @param df output from functions output_tp
#' @param path output folder location
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @return pymol script with residues colored based on average of uptake per residue.
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a<- output_tp(file_nm)
#' pymol_script_average_residue(df=a, replicates=3, pv_cutoff=0.01,
#' ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf), path=tempdir() )
#' pymol_script_average_residue(df=a,  path=tempdir())
#' }
#' @export
pymol_script_average_residue<-function(df,path="", ranges=c(-Inf, seq(-30, 30, by=10), Inf),
pv_cutoff=0.01, replicates=3){
  #####from HDX get data and

  oldwd<-getwd()
  on.exit(setwd(oldwd))
  setwd(path)

  for ( deut_time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(df[df$Deut.Time==deut_time,],  replicates)
    sd<-sd_timepoint(df[df$Deut.Time==deut_time,], replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut_time,], replicates)


    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(df1, pv, sd, pv_cutoff, replicates)

    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')

    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/(df1[,7])*10*cl1)###vector which has significant average

    ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df1[,end_col]), rep(0, length.out=max(df1[,end_col])*dim(fc.d)[2])))##make mock matrix
    #ac<-data.frame(matrix(ncol=dim(fc.d)[2] , nrow =max(df1[,end_col]), rep(0, length.out=dim(max(df1[,end_col]))*dim(fc.d)[2])))##make mock matrix
    for ( j in 1:dim(ac)[2]){
      for ( i in 1:dim(df1)[1]){

        sum.nm<-(ac[df1[i,start_col]:df1[i,end_col],j]+fc.d[i,j])
        ###create a vector of which is a sum of significance at position
        ac[df1[i,start_col]:df1[i,end_col],j]=sum.nm[1] ## assign new value to position
      }}


    ##prep of coverage
    coverage<-coverage_residue(df1, start_col, end_col)

    ave.p.cov<-ac/coverage ## sums of the significant avererages divided by coverage.
    ave.p.cov[ave.p.cov=="NaN"]<-0 ##remove NAN divisions introduced by division by no coverage


    #####preparation of average per residue data.frame, which will have

    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c("white", "grey45"))


    for ( i in 1:(length(xli)-1)){
      ave.p.cov[xli[i]< ave.p.cov & ave.p.cov < xli[i+1]] <- num_ass[i]}
    ave.p.cov[ave.p.cov==0]<- (-10000)

    si_apc<-abs(ave.p.cov)-9999
    cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
    si_apc<-si_apc*cv_mis+1

    ##assign name to colors in pallettes
    col_nm<-c("no_cov", "NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}

    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette

    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i],
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}
    ###write outputs per each state in the
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-9)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_average_", nm1[j],"_", deut_time, ".txt", sep="")
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0){
          res.txt<-c(res.txt, (paste(c("color ", col_nm[i],", resi ",
                                       pymol_str(which(si_apc[,j] ==i))), sep="", collapse="")))}}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color grey", "bg white",set_colors,res.txt ), fileConn)
      close(fileConn)}}

  leg_nm<-c("No coverage", "Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}

  pallette_ll(cbr1, leg_nm)
  return()}
