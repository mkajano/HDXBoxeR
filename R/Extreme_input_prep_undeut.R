
#' Makes input for Extreme for bimodal analysis.
#'
#' If data is missing it returns non-deuterated data in these columns.
#'
#' @param hm_dir directory in which all the folders which needs to be processed are
#' @param replicates number of replicates in sample
#' @param timepoints lists timepoints used in experiments.
#' @param output_path directory where output should be written
#' @return Inputs for extreme for all data prepared.
#' @examples
#' \donttest{
#' path_to_folders<-system.file("extdata",  package = "HDXBoxeR")

#' extreme_input_undeut(hm_dir=path_to_folders, replicates = 3,
#' timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#' }
#' @export
extreme_input_undeut<-function(hm_dir, replicates,timepoints, output_path="NA" ){

  oldwd<-getwd()
  on.exit(setwd(oldwd))

  if (output_path =="NA"){
    output_path<-hm_dir
  }

  setwd(hm_dir)
  dirs<-list.dirs(full.names = TRUE, recursive = TRUE)

  for (directories in dirs){
    message(paste0("processing:", directories))
    nd<-list.files(directories, pattern="Non-D")
    fd<-list.files(directories, pattern="Full-D")
    tp<-c(list.files(directories,pattern="s-"), list.files(directories,pattern="m-"), list.files(directories, pattern="h-"))
    if (all(is.na(c(nd, fd, tp)))==TRUE) {
      next         ## << Only line that differs
    } else {

      tp1<-as.data.frame(tp) %>% separate(tp, into = paste("V", 1:3, sep =""),  sep="([\\-])", extra="merge")


      hms=str_sub(tp1[,1], start=-1)
      nbs.t=as.numeric(str_sub(tp1[,1], end=-2))

      for ( i in 1:length(hms)){
        if (hms[i]=="h"){
          nbs.t[i]<-nbs.t[i]*3600
        } else if (hms[i]=="m"){
          nbs.t[i]<-nbs.t[i]*60
        } else (nbs.t[i]=nbs.t[i])
      }

      tp1[,1]<-nbs.t

      for (charge in unique(tp1[,3])){
        ind_one_charge=which(tp1[,3]==charge)
        one_rep<-tp1[which(tp1[,3]==charge),]
        one_nbs.t<-nbs.t[ind_one_charge]

        nd1<-as.data.frame(nd) %>% separate(nd, into = paste("V", 1:4, sep = c("-")), sep="([\\-])", extra="merge")
        fd1<-as.data.frame(fd) %>% separate(fd, into = paste("V", 1:4, sep = c("-")), sep="([\\-])", extra="merge")


        replicates_each<-as.data.frame(table(one_nbs.t))
        #sort_tp<-sort(tp1[,1])


        mock.df<-data.frame(c(NA, NA), c(NA, NA))
        colnames(mock.df)<-c("V1", "V2")
        ##preparation of files to read

        ###

        nm_nd<-c(paste("t-", str_sub(c(nd[which(nd1[,4]==charge)[1]]), end=-5), sep=""))
        nm_nd<-gsub("-", "_", nm_nd)
        if(nm_nd=="t_NA"){
          nm_nd<-"mock.df"
        } else {my.files <- c(nd[which(nd1[,4]==charge)[1]])
        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=FALSE)
          my.name <- nm_nd[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}
        }

        nm_fd<-c(paste("t-", str_sub(c(fd[which(fd1[,4]==charge)[1]]), end=-5), sep=""))
        nm_fd<-gsub("-", "_", nm_fd)
        if(nm_fd=="t_NA"){
          nm_fd<-"nd1"
        } else {my.files <- c(fd[which(fd1[,4]==charge)[1]])

        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=FALSE)
          my.name <- nm_fd[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}
        }




        nm1<-c(paste("t-", str_sub(c(tp), end=-5), sep="")[ind_one_charge] )
        nm1<-gsub("-", "_", nm1)
        my.files <- c( tp[ind_one_charge])
        for (k in 1:(length(my.files))){
          full_path=paste(directories,"/", my.files[k], sep="")
          # import the file
          cur.file <- read.csv(file = full_path, header=FALSE)
          my.name <- nm1[k]
          # assign the name to the object
          assign(paste(my.name), cur.file)}

        #### order the files correctly + add mock dif to lists of files that need to bound together.

        order_files<-orderv(one_rep)
        nm1_rb<-c()
        for ( ti in 1:length(timepoints)){
          nb_rep=length(which(one_rep[,1]==timepoints[ti]))
          nm_cb=which(one_rep[,1]==timepoints[ti])
          if (nb_rep == replicates ){
            nm1_rb<-c(nm1_rb, nm1[nm_cb])
          } else if (nb_rep<replicates & nb_rep >0){
            dif_mock=replicates-nb_rep
            nm1_rb<-c(nm1_rb, nm1[nm_cb])
            nm1_rb<-c(nm1_rb, rep(nm_nd, dif_mock))
          } else if (nb_rep==0){
            nm1_rb<-c(nm1_rb, rep(nm_nd, replicates))
          }
        }
        nm1_rb<-c(nm_nd,nm1_rb, nm_fd)

        dfB<-c()
        dfB<-mget(nm1_rb[1])[[1]]
        for ( rbi in 2:length(nm1_rb)){
          dfB<-qpcr.cbind.na(dfB, mget(nm1_rb[rbi])[[1]])}


        labs<-c(paste("undeut", sep=""))

        for ( ti in 1:length(timepoints)){
          for (tj in seq(replicates)){
            labs<-c(labs,paste(timepoints[ti],".0", tj, " sec", sep=""))}}
        labs<-c(labs,paste("TD",  sep=""))

        labsy<-rep(" ", length.out=length(labs))

        nms_df<-c()
        for ( ti in 1: length(labsy)){
          nms_df<-c(nms_df, labs[ti], labsy[ti])}

        colnames(dfB)<-nms_df

        dir_nm<-gsub("/", "", directories)

        output_name<-paste(output_path,"/", dir_nm,"_", charge,  ".csv", sep="")

        write.table(dfB, output_name,na = "",row.names = FALSE, sep = ",")}}

  }}


