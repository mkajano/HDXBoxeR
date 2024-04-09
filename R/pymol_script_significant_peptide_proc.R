#' Writes a text files with pymol scripts to list significant peptides
#'
#' Function write a script that can be used in pymol to color structure.
#' Number of colors and corresponding to them ranges can be defined by user.
#'
#' @param input_proc Dataframe with organized procent deuteration data. Input generated using output_tp(, percent=T) function.
#' @param input_up Dataframe with organized deuteration uptake. Input generated using output_tp() function.
#' @param path location where the Pymol scripts will be saved
#' @param pv_cutoff p-value cutoff here set up to 0.01
#' @param replicates number of replicates in sample. Default set to 3.
#' @param ranges ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30, by=10), Inf)
#' @param order.pep flag allowing to either order peptide acccording to the peptide length (default), or to position in the protein sequence.
#' @return pymol script with colors assigned per peptide
#' @examples
#' \donttest{
#' file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a_up<- output_tp(file_nm)
#' a_proc<- output_tp(file_nm, percent=TRUE)
#' pymol_script_significant_peptide_proc(input_proc=a_proc,
#' input_up=a_up,  path=tempdir(),replicates=3, pv_cutoff=0.01,
#' ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf), order.pep=TRUE)
#' }
#' @export
pymol_script_significant_peptide_proc<-function(input_proc,input_up, path="", ranges=c(-Inf, seq(-30, 30, by=10), Inf),
                                                pv_cutoff=0.01, replicates=3, order.pep=TRUE){
  dfup=input_up
  df=input_proc

  oldwd<-getwd()
  on.exit(setwd(oldwd))
  setwd(path)

  #####from HDX get data and
  for ( deut.time in(unique(df$Deut.Time))){
    pv<-pv_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    sd<-sd_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    dfu<-ave_timepoint(dfup[dfup$Deut.Time==deut.time,], replicates)
    df1<-ave_timepoint(df[df$Deut.Time==deut.time,], replicates)

    #preparation significance per residue & coverage
    cl1<-significant_peptide_uptake(dfu, pv, sd, pv_cutoff, replicates)
    start_col<-which(colnames(df1)=='Start')
    end_col<-which(colnames(df1)=='End')

    ###preparation of the coloring ranges. Difference of average.
    fc.d<-data.frame((df1[,7]-df1[,8:dim(df1)[2]])/10*cl1)###vector which has significant average
    ##sumarized occurances of peptides
    si.f=fc.d
    #####preparation of average per residue data.frame, which will have
    xli=ranges/10; num_ass<-c(-10001:(-10001-(length(xli)-2)))
    for ( i in 1:(length(xli)-1)){
      si.f[xli[i]<  si.f &  si.f < xli[i+1]] <- num_ass[i]}
    si.f[ si.f==0]<- (-10000)
    si_apc<-abs(si.f)-9999



    cbr1<-color_ranges_Blue_Red_heat_map(ranges=xli, c( "white"))

    ##assign name to colors in pallettes
    col_nm<-c("NSig")
    for ( i in 1:(length(ranges)-1)){
      col_nm<-c(col_nm, paste("col_", ranges[i],"_", ranges[i+1], sep=""))}


    rgb_col<-col2rgb(cbr1, alpha = FALSE) ## function to return rgb values for colors in pallette

    ##set_color 0%, [0 , 0 , 120], make command for pymol, to set colors
    set_colors<-c()
    for ( i in 1:length(col_nm)){
      set_colors<-c(set_colors, paste("set_color ", col_nm[i],", [",  rgb_col[1,i],
                                      ",", rgb_col[2,i], ",", rgb_col[3,i], "]", sep=""))}


    ###write outputs per each state in the
    nm1<-str_sub(colnames(df1[8:dim(df1)[2]]), start=4, end=-8)
    for (j in 1:dim(si_apc)[2]){
      output_name<-paste("pymol_all_peptides_proc_", nm1[j],"_", deut.time, ".txt", sep="")

      res.txt<-c()
      len.pep<-c()
      for ( i in 1:length(col_nm)){
        if (length(which(si_apc[,j]==i)) !=0 ) {
          pep_nb<-which(si_apc[,j] == i)
          for ( k in pep_nb){line<-c()
          line<- paste(c("color ", col_nm[i],", resi ", df1[k,start_col], "-", df1[k,end_col] , sep=""))
          len.pep<-c(len.pep, df1[k,end_col]-df1[k,start_col] )
          res.txt<-c(res.txt,
                     paste(line, sep="' '", collapse=""))}
        }}

      if (order.pep==T){
        res.txt<-res.txt[rev(order(len.pep))]
        message("peptides ordered according to peptide length")
      } else if (order.pep==FALSE){
        message("peptides ordered according to position in sequence")}

      fileConn<-file(output_name)
      writeLines(c("hide","show cartoon","color black", "bg white",
                   set_colors, res.txt ), fileConn)
      close(fileConn)}}
  leg_nm<-c("Not Sig")
  for ( i in 1:(length(ranges)-1)){
    leg_nm<-c(leg_nm, paste(ranges[i],":", ranges[i+1], "%", sep=""))}

  pallette_ll(cbr1, leg_nm)
  return()}
