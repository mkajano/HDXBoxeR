#' Writes a text files with pymol scripts to list significant residues.
#'
#' Function write a script that can be used in pymol to color structure.
#' Number of colors and corresponding to them ranges can be defined by user.
#' Residues are colored by average procent_deuteration from the significant peptides per residues.
#'
#' @details Percent-deuteration counterpart of [pymol_script_significant_residue()]: significance is
#'   tested on the uptake table `input_up` (Da, the package-wide rule), while
#'   residues are colored by the largest significant difference in percent deuteration, first state minus each state (from `input_proc`); [pymol_script_significant_residue()] uses the relative change in uptake.
#'
#' @param input_up Dataframe with organized deuteration uptake. Input generated using output_tp() function.
#'   Note: in this function `input_up` is the first argument (as in earlier versions); name both arguments to be safe.
#' @param input_proc Dataframe with organized procent deuteration data. Input generated using output_tp(percent = TRUE).
#' @param path location where the Pymol scripts will be saved. Default `""`: the working directory.
#' @param alpha cutoff for critical interval
#' @param pvalue p-value cutoff here set up to 0.01
#' @param pv_cutoff deprecated, use `pvalue`.
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param p_adjust multiple-testing correction of the Welch p-values, passed to [stats::p.adjust()]: `"none"` (default, raw p-values) or e.g. `"BH"`. A peptide is significant when its (adjusted) p-value is below the cutoff AND its difference exceeds the critical interval.
#' @param resno_offset integer added to every residue number written to the PyMOL script, so the data's numbering matches the structure's (e.g. `-20` when the structure numbering starts 20 residues later, as with an N-terminal tag). Default `0`.
#' @return pymol script with residues colored based on average of procent deuteration per residue.
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a_up<- output_tp(file_nm)
#' a_proc<- output_tp(file_nm, percent=TRUE)
#' pymol_script_significant_residue_proc(input_proc=a_proc,
#' input_up=a_up, path=tempdir(), replicates=3, alpha=0.01, pvalue=0.01,
#' ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf))
#' }
#' @export
pymol_script_significant_residue_proc<-function(input_up,input_proc,path="",
                                                ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                                                alpha=0.01,
                                                pvalue=0.01, replicates=3, p_adjust="none", resno_offset=0, pv_cutoff){
  if (!missing(pv_cutoff)) { .hdx_deprecated("pv_cutoff", "pvalue"); pvalue <- pv_cutoff; if (missing(alpha)) alpha <- pv_cutoff }
  ## input_up comes first here (as in 0.0.2), unlike the other _proc
  ## functions; tables given the other way round are swapped back
  is_pct <- function(d) any(grepl("_Deut\\.\\._[0-9]+$", colnames(d)))
  if (is_pct(input_up) && !is_pct(input_proc)) {
    warning("input_up looks like a percent table and input_proc like an uptake table: ",
            "they were swapped. Name the arguments to avoid this.", call. = FALSE)
    tmp <- input_up; input_up <- input_proc; input_proc <- tmp
  }
  .hdx_check_aligned(input_proc, input_up)

  oldwd<-getwd()
  on.exit(setwd(oldwd))
  if (nzchar(path)) setwd(path)

   #####from HDX get data and
  dfup=input_up
  df=input_proc
  for ( deut.time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = replicates, p_adjust = p_adjust)
    sd<-sd_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = replicates)
    dfu<-ave_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates = replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut.time,], replicates = replicates)

    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(dfu, pv, sd, alpha, pvalue,replicates)

    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')

    ###preparation of the coloring ranges. Difference of average.
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/10*cl1)###vector which has significant average
    ###per residue maximum uptake value is chosen. ->
    max.ac1<-c()
    for ( j in 1:dim(fc.d)[2]){
      ac<-c()
      ac1<-c()
      ac2<-c()
      for ( i in 1:dim(df1)[1]){
        ac<-rep(0, length=max(df1[,end_col]))
        ac[df1[i,start_col]:df1[i,end_col]]=fc.d[i,j]
        ##make multiple vectors which have 1 at position which peptide covers
        ac1<-c(ac1, ac)}
      ac2=data.frame(matrix(ac1, nrow =dim(df1)[1], byrow=T))
      max.a<-c()
      for ( k in 1:dim(ac2)[2]){
        ind1<-which.max(abs(ac2[,k]))
        nb1<-(ac2[ind1,k])
        max.a<-c(max.a, nb1)}
      max.ac1<-c(max.ac1, max.a)}
    max.ac2=data.frame(matrix(max.ac1, ncol = dim(fc.d)[2]))


    ##prep of coverage
    coverage<-coverage_residue(df1, start_col, end_col)

    #####preparation of average per residue data.frame, which will have
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      max.ac2[xli[i]< max.ac2 & max.ac2 < xli[i+1]] <- num_ass[i]}
    max.ac2[max.ac2==0]<- (-10000)

    si_apc<-abs(max.ac2)-9999
    cv_mis=coverage; cv_mis[cv_mis > 1]<- (1)###define lack of coverage
    si_apc<-si_apc*cv_mis+1

    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c( "grey45", "white"))


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
    nm1<-sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", sub("^av_", "", colnames(df1)[8:dim(df1)[2]]))
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_maxuptake_proc_", nm1[j],"_", deut.time, ".txt", sep="")
      res.txt<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0){
          res.txt<-c(res.txt, (paste(c("color ", col_nm[i],", resi ",
                                       pymol_str(which(si_apc[,j] ==i) + resno_offset)), sep="", collapse="")))}}
      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color grey", "bg white",set_colors,res.txt ), fileConn)
      close(fileConn)}}

  leg_nm<-c("No coverage", "Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}

  pallette_ll(cbr1, leg_nm)
  return()}
