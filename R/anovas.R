# ## work in progress +tests.
# library(HDXBoxeR)
# library(RColorBrewer)
# library(stringi)
# library(stringr)
# library(ggplot2)
# library(car)
# library(dplyr)
# library(multcomp)
# library(tidyr)
#
# library(ggstatsplot)
#
# #### variables
# path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/all_data.csv"
#
# names_states<- nm_states(path)
# set.seed(123)
#
# end_ntr<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)
# end_acd<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)+80
#
#
#
# ####functions
# boxer_input_to_long<-function(path, states, seq.match=F,...){
#
#   a5b<-output_tp(path, states =states, ...)
#
#   v5 <- gather(a5b, Experiment, X..Deut,
#                names(a5b)[7]:names(a5b)[length(names(a5b))], factor_key=TRUE)
#   nms<-data.frame(str_split_fixed(v5$Experiment, 'X', 2))
#   v5$Protein.State <- str_sub(nms[,1], end = -2)
#   v5$Experiment<-paste("X", nms[,2], sep="")
#
#   v5$nb_prot<-nb_exch_deut(v5)
#   v5$frac_deut<-as.numeric(v5$X..Deut)/v5$nb_prot
#
#   return(v5)
# }
#
# plt_gg_btwstats<-function(df, names){
#   inds<-c()
#   for (i in 1: length(names)){
#     inds<-c(inds, which(df$Protein.State==names[i]))}
#
#   ggbetweenstats(
#     data = df[inds,],
#     x = Protein.State,
#     y = frac_deut,
#     type = "nonparametric", # ANOVA or Kruskal-Wallis
#     var.equal = FALSE, # ANOVA or Welch ANOVA
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = TRUE,
#     bf.message = TRUE,
#     p.adjust.method = "holm" ,
#     package      = "ggsci",
#     palette      = "default_igv"
#   )
# }
#
#
#
# plt_gg_groupedbtw<-function(df, names){
#   inds<-c()
#   for (i in 1: length(names)){
#     inds<-c(inds, which(df$Protein.State==names[i]))}
#
#   grouped_ggbetweenstats(
#     data = df[inds,],
#     x = Protein.State,
#     y = frac_deut,
#     grouping.var = domain,
#     type = "nonparametric", # ANOVA or Kruskal-Wallis
#     var.equal = FALSE, # ANOVA or Welch ANOVA
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = TRUE,
#     bf.message = TRUE,
#     p.adjust.method = "holm" ,
#     package      = "ggsci",
#     palette      = "default_igv"
#   )
# }
#
# plt_gg_groupedwithin<-function(df, names){
#   inds<-c()
#   for (i in 1: length(names)){
#     inds<-c(inds, which(df$Protein.State==names[i]))}
#
#   grouped_ggwithinstats(
#     data = df[inds,],
#     x = Protein.State,
#     y = frac_deut,
#     grouping.var     = domain,
#     type = "nonparametric", # ANOVA or Kruskal-Wallis
#     var.equal = FALSE, # ANOVA or Welch ANOVA
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = TRUE,
#     bf.message = TRUE,
#     p.adjust.method = "holm" ,
#     package      = "ggsci",
#     palette      = "default_igv"
#   )
# }
#
# plt_gg_withinstats<-function(df, names){
#   inds<-c()
#   for (i in 1: length(names)){
#     inds<-c(inds, which(df$Protein.State==names[i]))}
#
#   ggwithinstats(
#     data = df[inds,],
#     x = Protein.State,
#     y = frac_deut,
#     type = "nonparametric", # ANOVA or Kruskal-Wallis
#     var.equal = FALSE, # ANOVA or Welch ANOVA
#     plot.type = "box",
#     pairwise.comparisons = TRUE,
#     pairwise.display = "significant",
#     centrality.plotting = TRUE,
#     bf.message = TRUE,
#     p.adjust.method = "holm" ,
#     package      = "ggsci",
#     palette      = "default_igv"
#   )
# }
#
# added_3domains<-function(df, end_ntr, end_acd){
#
#
#   bt<-df
#   uq_st<-unique(bt$Protein.State)
#   bt$domain<-"CTR"
#   for (st in 1:length(uq_st)){
#     ind_to_be_changed<-c()
#     ind_to_be_changed<-which(bt$Protein.State==uq_st[st] & bt$Start<=end_acd[st]& bt$End<=end_acd[st])
#     bt$domain[ind_to_be_changed]<-"ACD"
#     ind_to_be_changed<-c()
#     ind_to_be_changed<-which(bt$Protein.State==uq_st[st] & bt$End<=end_ntr[st])
#     bt$domain[ind_to_be_changed]<-"NTR"
#   }
#
#
#   indn<-which(bt$domain=="NTR")
#   inda<-which(bt$domain=="ACD")
#   indc<-which(bt$domain=="CTR")
#   bn<-bt[indn,]
#   ba<-bt[inda,]
#   bc<-bt[indc,]
#   return(list(bt, bn, ba, bc))
# }
# all_data_long_input<-function(path){
#   b1<-read.csv(path)
#   bt<-b1[which(b1$Deut.Time != "0s" & b1$Deut.Time != "FD" ),]
#   bt$nb_prot<-nb_exch_deut(bt)
#   bt$frac_deut<-as.numeric(bt$X..Deut)/bt$nb_prot
#   return(bt)
# }
#
#
# ####
# path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/all_data.csv"
#
# set.seed(123)
#
# end_ntr<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)
# end_acd<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)+80
#
#
#
# uq_st<-nm_states(path)
# b<-all_data_long_input(path)
# bt<-added_3domains(b, end_ntr, end_acd)[[1]]
# bn<-added_3domains(b, end_ntr, end_acd)[[2]]
# ba<-added_3domains(b, end_ntr, end_acd)[[3]]
# bc<-added_3domains(b, end_ntr, end_acd)[[4]]
#
# ##for groupped analysis //control_vs_bump
# l5<-boxer_input_to_long(path, uq_st[5:6])
# v5<-added_3domains(l5, end_ntr, end_acd)[[1]]
#
# l1<-boxer_input_to_long(path, uq_st[c(1,4)])
# v1<-added_3domains(l1, end_ntr, end_acd)[[1]]
#
# l151<-boxer_input_to_long(path, uq_st[7:8])
# v151<-added_3domains(l151, end_ntr, end_acd)[[1]]
#
# l515<-boxer_input_to_long(path, uq_st[9:10])
# v515<-added_3domains(l515, end_ntr, end_acd)[[1]]
#
# #### Chimeras
#
#
# p1<-plt_gg_btwstats(bt, uq_st[c(2,5,7,9)])
# p2<-plt_gg_btwstats(bn, uq_st[c(2,5,7,9)])
# p3<-plt_gg_btwstats(ba, uq_st[c(2,5,7,9)])
# p4<-plt_gg_btwstats(bc, uq_st[c(2,5,7,9)])
#
#
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 2),
#   annotation.args = list(
#     title = "Comparison of the fraction deutered per protein fragment"
#   )
# )
#
# ####bumps not paired. all
#
# p1<-plt_gg_btwstats(bt, uq_st[c(1,4)])
# p2<-plt_gg_btwstats(bt, uq_st[c(5,6)])
# p3<-plt_gg_btwstats(bt, uq_st[c(7,8)])
# p4<-plt_gg_btwstats(bt, uq_st[c(9,10)])
#
#
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 1),
#   annotation.args = list(
#     title = "Comparison of all peptides, wild-type to bump"
#   )
# )
#
# ####
#
#
# ####bumps not paired NTR
#
# p1<-plt_gg_btwstats(bn, uq_st[c(1,4)])
# p2<-plt_gg_btwstats(bn, uq_st[c(5,6)])
# p3<-plt_gg_btwstats(bn, uq_st[c(7,8)])
# p4<-plt_gg_btwstats(bn, uq_st[c(9,10)])
#
#
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 1),
#   annotation.args = list(
#     title = "Comparison of NTR peptides, wild-type to bump"
#   )
# )
#
#
# ###
#
# ####bumps not paired. ACD
#
# p1<-plt_gg_btwstats(ba, uq_st[c(1,4)])
# p2<-plt_gg_btwstats(ba, uq_st[c(5,6)])
# p3<-plt_gg_btwstats(ba, uq_st[c(7,8)])
# p4<-plt_gg_btwstats(ba, uq_st[c(9,10)])
#
#
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 1),
#   annotation.args = list(
#     title = "Comparison of ACD peptides, wild-type to bump"
#   )
# )
#
# ### paired all
#
# p1<-plt_gg_withinstats(v1, unique(v1$Protein.State))
# p2<-plt_gg_withinstats(v5, unique(v5$Protein.State))
# p3<-plt_gg_withinstats(v151, unique(v151$Protein.State))
# p4<-plt_gg_withinstats(v515, unique(v515$Protein.State))
#
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 1),
#   annotation.args = list(
#     title = "Comparison of all peptides, wild-type to bump, paired"
#   )
# )
#
# p1<-plt_gg_groupedwithin(v1, unique(v1$Protein.State))
#
# p2<-plt_gg_groupedwithin(v5, unique(v5$Protein.State))
# p3<-plt_gg_groupedwithin(v151, unique(v151$Protein.State))
# p4<-plt_gg_groupedwithin(v515, unique(v515$Protein.State))
#
#
# combine_plots(
#   list(p1, p2),
#   plotgrid.args = list(nrow = 2),
#   annotation.args = list(
#     title = "Comparison of all peptides, wild-type to bump, paired"
#   )
# )
#
#
# combine_plots(
#   list(p3, p4),
#   plotgrid.args = list(nrow = 2),
#   annotation.args = list(
#     title = "Comparison of all peptides, wild-type to bump, paired"
#   )
# )
