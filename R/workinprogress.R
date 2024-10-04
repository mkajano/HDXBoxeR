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
# sd_ci_timepoint<-function(df, replicates=3, alpha=0.01) {
#   ### calculate standard deviation of sample 1 for all data points
#   nb_sets=(dim(df)[2]-6)/replicates
#   nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
#   nm_root<-str_sub(nm_root, end = -3)
#   sd_nm<-paste("sd_ci_", nm_root, sep="")
#
#   sd1<-c(); for ( j in 1:nb_sets) {
#     for (i in 1:dim(df)[1]) {x1<-df[i,7:(7+replicates-1)];
#     x2<-df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]
#
#     tvalue=abs(qt(alpha, replicates*2-2))
#     tt<-t.test(x1, x2, conf.level = c(1-alpha))
#     sd.ci<-(tt$conf.int[2]-tt$conf.int[1])/tvalue*sqrt(replicates)
#
#
#     sd1<-c(sd1, sd.ci)}
#   }
#
#   sd2<-data.frame(matrix(sd1, ncol=nb_sets , byrow = FALSE))
#   colnames(sd2)<-sd_nm
#   sd2<-data.frame(df[,1:6], sd2)
#   return(sd2)}
#
#
# path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/all_data.csv"
#
#
# path21="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/heterooligomers/B1_heterooligomers_Allresults_culled.csv"
# path25="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/heterooligomers/B5_heterooligomers_Allresults_culled.csv"
#
#
#
# names_states<- nm_states(path)
# a5<-output_tp(path, states =names_states[c(5:6)], seq_match = F )
# a1<-output_tp(path, states =names_states[c(1,4)], seq_match = F )
#
# a515<-output_tp(path, states =names_states[c(9:10)], seq_match = F )
# a151<-output_tp(path, states =names_states[c(7:8)], seq_match = F )
#
# a5a<-output_tp(path, states =names_states[c(5)], seq_match = F )
# a1a<-output_tp(path, states =names_states[c(2)], seq_match = F )
#
# a515a<-output_tp(path, states =names_states[c(9)], seq_match = F )
# a151a<-output_tp(path, states =names_states[c(7)], seq_match = F )
#
#
#
#
# names_states21<- nm_states(path21)
# h1<-output_tp(path21, states =names_states21[c(1:2)], seq_match = F, times = "5.00s")
#
# names_states25<- nm_states(path25)
# h5<-output_tp(path25, states =names_states25[c(1:2)], seq_match = F, times = "5.00s" )
#
#
#
#
# nb_ex1<-nb_exch_deut(a1)
# av1<-ave_timepoint(a1)
# sd1<-sd_timepoint(a1)
# p1<-pv_timepoint(a1)
# dav<-dif_ave(av1)
# sd_ci<-sd_ci_timepoint(a1, alpha=0.05)
# d1<-dim(av1)[2]
#
#
# boxplot(c(a1[,7:9],a5[,7:9]))
#
# #SD=sqrt(replicates)*(uper_limit-lower_limit)/t
# # https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals5.html
# # https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
#
# replicates=3
# alpha=0.01
#
# sd_ci<-sd_ci_timepoint(a1)
#
# tvalue=abs(qt(alpha, replicates*2-2))
# ts<-t.test(a1[15,7:9], a1[15,10:12], conf.level = c(1-alpha))
# sd.ci<-(ts$conf.int[2]-ts$conf.int[1])/tvalue*sqrt(replicates)
#
# t5<-output_tc(path, states =names_states[c(5:6)], seq_match = F )
# t1<-output_tc(path, states =names_states[c(1,4)], seq_match = F )
#
#
# woods_CI_plot_frac(t5,replicates = 3)
# deuteration_woods_timepoints_frac(a1)
# deuteration_woods_timepoints_frac(a5)
#
#
# deuteration_woods_timepoints_frac(a1, replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac(a5, replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac(a151, replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac(a515, replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac_list(list(a1a, a151a), replicates = 3, ylim = c(0,0.8))
# deuteration_woods_timepoints_frac_list(list(a5a, a515a), replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac_list(list(a1a, a151a, a5a), replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac_list(list(a1a, a151a, h1), replicates = 3, ylim = c(0,0.8))
#
# deuteration_woods_timepoints_frac_list(list(a5a, a515a, h5), replicates = 3, ylim = c(0,0.8))
#
#
# deuteration_woods_timepoints_frac_list(list(h5), replicates = 3, ylim = c(0,0.8))
# deuteration_woods_timepoints_frac_list(list(h1), replicates = 3, ylim = c(0,0.9))
#
# deuteration_woods_timepoints_frac_list(list(a1a), replicates = 3, ylim = c(0,0.8))
# deuteration_woods_timepoints_frac_list(list(a151a), replicates = 3, ylim = c(0,0.8))
#
#
# deuteration_woods_timepoints_frac_list(list(a5a), replicates = 3, ylim = c(0,0.8))
# deuteration_woods_timepoints_frac_list(list(a515a), replicates = 3, ylim = c(0,0.8))
#
# pymol_script_significant_peptide(a1, "C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/")
#
# pymol_script_significant_peptide(a5, "C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/")
#
#
# plot_peptide_sig_tp(a5)
# plots_vol_tp(a1)
#
# #### anova tests
#
#
# uq_st<-unique(bt$Protein.State)
# end_ntr<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)
#
# end_ntr<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)
# end_acd<-c(85, 85, 85, 85, 75, 75, 85,85, 75, 75)+80
#
# ###2 domains taken into account
#
# bt$domain<-"ACD_CTR"
# for (st in 1:length(uq_st)){
#   ind_to_be_changed<-c()
#   ind_to_be_changed<-which(bt$Protein.State==uq_st[st] & bt$End<=end_ntr[st])
#   bt$domain[ind_to_be_changed]<-"NTR"
# }
#
# indn<-which(bt$domain=="NTR")
# indac<-which(bt$domain=="ACD_CTR")
# bn<-bt[indn,]
# bac<-bt[indac,]
#
# indn<-which(bt$domain=="NTR")
# indac<-which(bt$domain=="ACD_CTR")
# bn<-bt[indn,]
# bac<-bt[indac,]
#
# p <- ggbetweenstats(data = bt[states1,],
#                     x = Protein.State,
#                     y = frac_deut )
# p
#
# # extracting details from statistical tests
# extract_stats(p)
#
#
# ## 3 domains
#
#
# bt$domain<-"CTR"
# for (st in 1:length(uq_st)){
#   ind_to_be_changed<-c()
#   ind_to_be_changed<-which(bt$Protein.State==uq_st[st] & bt$Start<=end_acd[st]& bt$End<=end_acd[st])
#   bt$domain[ind_to_be_changed]<-"ACD"
#   ind_to_be_changed<-c()
#   ind_to_be_changed<-which(bt$Protein.State==uq_st[st] & bt$End<=end_ntr[st])
#   bt$domain[ind_to_be_changed]<-"NTR"
# }
#
#
# indn<-which(bt$domain=="NTR")
# inda<-which(bt$domain=="ACD")
# indc<-which(bt$domain=="CTR")
# bn<-bt[indn,]
# ba<-bt[inda,]
# bc<-bt[indc,]
#
#
#
#
#
# all_v5<-added_3domains(v5, c(65, 65), c(155, 155))
#
# vb<-all_v5[[1]]
# plt_gg_groupedwithin(vb, uq_st[c(5,6)])
# plt_gg_withinstats(vb, uq_st[c(5,6)])
#
# plt_gg_btwstats(vb, uq_st[c(5,6)])
# plt_gg_groupedbtw(bt, uq_st[c(5,6)])
# ####
#
# # bumps plots t-test
#
# combine_plots(
#   list(p1, p2, p3, p4),
#   plotgrid.args = list(nrow = 2),
#   annotation.args = list(
#     title = "Comparison of life expectancy between 1957 and 2007",
#     caption = "Source: Gapminder Foundation"
#   )
# )
#
#
# plt_gg_groupedstats(bt, uq_st[c(4,1)])
#
# plt_gg_groupedstats(v5, uq_st[c(5,6)])
#
# plt_gg_btwstats(bt, uq_st[c(4,1)])
# plt_gg_btwstats(bn, uq_st[c(4,1)])
# plt_gg_btwstats(ba, uq_st[c(4,1)])
#
#
# plt_gg_btwstats(bt, uq_st[c(5,6)])
# plt_gg_btwstats(bn, uq_st[c(5,6)])
# plt_gg_btwstats(ba, uq_st[c(5,6)])
#
#
#
# plt_gg_btwstats(bt, uq_st[c(1,2,4)])
# plt_gg_btwstats(bn, uq_st[c(1,2,4)])
# plt_gg_btwstats(ba, uq_st[c(1,2,4)])
#
# plt_gg_btwstats(bt, uq_st[c(1,2,3,4)])
# plt_gg_btwstats(bn, uq_st[c(1,2,3,4)])
#
#
# plt_gg_btwstats(bt, uq_st[c(2,3,4)])
# plt_gg_btwstats(bn, uq_st[c(1,2,4)])
# plt_gg_btwstats(ba, uq_st[c(1,2,4)])
#
#
# plt_gg_btwstats(bt, uq_st[c(1,4)])
# plt_gg_btwstats(bn, uq_st[c(1,4)])
# plt_gg_btwstats(ba, uq_st[c(1,4)])
#
# plt_gg_btwstats(bt, uq_st[c(2,7,8,5,6,9,10)])
#
# plt_gg_btwstats(bn, uq_st[c(2,7,8,5,6,9,10)])
# plt_gg_btwstats(bac, uq_st[c(2,7,8,5,6,9,10)])
# plt_gg_btwstats(ba, uq_st[c(2,7,8,5,6,9,10)])
# plt_gg_btwstats(bc, uq_st[c(2,7,8,5,6,9,10)])
#
#
#
# plt_gg_btwstats(bt, uq_st[c(5,6,9,10)])
#
# plt_gg_btwstats(bn, uq_st[c(5,6,9,10)])
# plt_gg_btwstats(ba, uq_st[c(5,6,9,10)])
# plt_gg_btwstats(bac, uq_st[c(5,6,9,10)])
#
#
#
# plt_gg_btwstats(bt, uq_st[c(2,7,8)])
#
# plt_gg_btwstats(bn, uq_st[c(2,7,8)])
# plt_gg_btwstats(ba, uq_st[c(2,7,8)])
# plt_gg_btwstats(bac, uq_st[c(2,7,8)])
#
#
# plt_gg_btwstats(ba, uq_st[c(2,7,9)])
# plt_gg_btwstats(ba, uq_st[c(5,7,9)])
# plt_gg_btwstats(ba, uq_st[c(2,5,7,9)])
#
#
#
# plt_gg_btwstats(bn, uq_st[c(2,7,9)])
# plt_gg_btwstats(bn, uq_st[c(5,7,9)])
# plt_gg_btwstats(bn, uq_st[c(2,5,7,9)])
#
#
#
# plt_gg_btwstats(bt, uq_st[c(2,5,7,9)])
#
# plt_gg_btwstats(bt, uq_st[c(2,7,9)])
# plt_gg_btwstats(bt, uq_st[c(5,7,9)])
#
#
#
#
# states1<-which(bn$Protein.State==uq_st[c(2,7,8)])
# states1c<-which(bn$Protein.State==uq_st[c(7,8)])
# states1b<-which(bn$Protein.State==uq_st[c(1,4)])
# states5<-which(bn$Protein.State==uq_st[c(5,6,9,10)])
#
#
# states15<-which(bn$Protein.State==uq_st[c(2,7,8,5,6,9,10)])
#
#
# ggbetweenstats(
#   data = bt,
#   x = Protein.State,
#   y = frac_deut,
#   type = "parametric", # ANOVA or Kruskal-Wallis
#   var.equal = TRUE, # ANOVA or Welch ANOVA
#   plot.type = "box",
#   pairwise.comparisons = TRUE,
#   pairwise.display = "significant",
#   centrality.plotting = FALSE,
#   bf.message = FALSE
# )
#
#
#
#
#
#
#
#
#
# group_by(bt, Protein.State) %>%
#   summarise(
#     mean = mean(frac_deut, na.rm = TRUE),
#     sd = sd(frac_deut, na.rm = TRUE)
#   )
#
#
# oneway.test(frac_deut ~ Protein.State,
#             data = bn,
#             var.equal = FALSE # assuming unequal variances
# )
#
#
# res_aov <- aov(frac_deut ~ Protein.State,
#                data = bn
# )
#
# summary(res_aov)
#
# ###
#
# ggbetweenstats(
#   data = bn[states1c,],
#   x = Protein.State,
#   y = frac_deut,
#   type = "parametric", # ANOVA or Kruskal-Wallis
#   var.equal = TRUE, # ANOVA or Welch ANOVA
#   plot.type = "box",
#   pairwise.comparisons = TRUE,
#   pairwise.display = "significant",
#   centrality.plotting = TRUE,
#   bf.message = TRUE
# )
#
# library(multcomp)
#
# # Tukey HSD test:
# post_test <- glht(res_aov,
#                   linfct = mcp(Protein.State = "Tukey")
# )
#
# summary(post_test)
#
# ggplot(bt[indn,]) +
#   aes(x = Experiment, y = frac_deut, color = Experiment) +
#   geom_jitter() +
#   theme(legend.position = "none")
#
# boxplot(frac_deut ~ Protein.State,
#         data = bt[indn,]
#         )
#
# kruskal.test(frac_deut ~ Experiment,
#              data = bt
# )
# leveneTest(frac_deut ~ Protein.State,
#            data = bt
# )
#
# res_aov <- aov(frac_deut ~ Experiment,
#                data = bt
# )
#
# summary(res_aov)
#
# par(mfrow = c(1, 2)) # combine plots
#
# # 1. Homogeneity of variances
# plot(res_aov, which = 3)
#
# # 2. Normality
# plot(res_aov, which = 2)
#
# par(mfrow = c(1, 2)) # combine plots
#
# # histogram
# hist(res_aov$residuals)
#
# # QQ-plot
# library(car)
# qqPlot(res_aov$residuals,
#        id = FALSE # id = FALSE to remove point identification
# )
#
# shapiro.test(res_aov$residuals)
#
#
#
# ggplot(bt) +
#   aes(x = Experiment, y = frac_deut) +
#   geom_boxplot()
#
#
#
#
#
# #
# plots_diff_tp(a1)
#
# plot(av1[,7]/nb_ex1, pch=20)
# points(av1[,10]/nb_ex1, col=5, pch=20)
#
#
# plot(av1[,9]/nb_ex1, pch=20)
# points(av1[,10]/nb_ex1, col=5, pch=20)
#
#
# plot(av1[,8]/nb_ex1, col=2, pch=20)
# points(av1[,9]/nb_ex1, col=4, pch=20)
#
# plot(av1[,7]/nb_ex1, col=2, pch=20)
# points(av1[,8]/nb_ex1, col=4, pch=20)
#
#
# plot(av1[,7]/nb_ex1, av1[,8]/nb_ex1)
#
#
# plot(av1[,7]/nb_ex1, av1[,10]/nb_ex1, pch=20)
#
# plot(av1[,7]/nb_ex1, av1[,10]/nb_ex1, pch=20, xlab="", ylab="")
# points(c(0,1),      c(0,1), lty=2, type="l")
# points(c(0,1), c(0,1)+0.04, lty=2, type="l")
# points(c(0,1), c(0,1)-0.04, lty=2, type="l")
#
#
#
#
# head(av1) ###average differences (against first Protein State in the file)
# da1<-dif_ave(av1)
#
# CI_single(sd1[,7], 3)
# CI_2pts(sd1[,7], sd1[,10],3)
#
#
# nb_ex5<-nb_exch_deut(a5)
#
#
# av5<-ave_timepoint(a5) sd5<-sd_timepoint(a5)
#
# plot(av5[,7]/nb_ex5, av5[,8]/nb_ex5, pch=20, xlab="", ylab="") points(c(0,1),
#                                                                       c(0,1), lty=2, type="l") points(c(0,1), c(0,1)+0.04, lty=2, type="l")
# points(c(0,1), c(0,1)-0.04, lty=2, type="l")
#
# plot(av5[,7]/nb_ex5, pch=20) points(av5[,8]/nb_ex5, col=2, pch=20)
# points(av5[,9]/nb_ex5, col=4, pch=20) points(av5[,10]/nb_ex5, col=5, pch=20)
#
