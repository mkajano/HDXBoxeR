## work in progress +tests.
library(HDXBoxeR)
library(RColorBrewer)
library(stringi)
library(stringr)
library(ggplot2)
library(car)


sd_ci_timepoint<-function(df, replicates=3, alpha=0.01) {
  ### calculate standard deviation of sample 1 for all data points
  nb_sets=(dim(df)[2]-6)/replicates
  nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
  nm_root<-str_sub(nm_root, end = -3)
  sd_nm<-paste("sd_ci_", nm_root, sep="")

  sd1<-c(); for ( j in 1:nb_sets) {
    for (i in 1:dim(df)[1]) {x1<-df[i,7:(7+replicates-1)];
    x2<-df[i,(6+(j-1)*replicates+1):(6+(j-1)*replicates+replicates)]

    tvalue=abs(qt(alpha, replicates*2-2))
    tt<-t.test(x1, x2, conf.level = c(1-alpha))
    sd.ci<-(tt$conf.int[2]-tt$conf.int[1])/tvalue*sqrt(replicates)


    sd1<-c(sd1, sd.ci)}
  }

  sd2<-data.frame(matrix(sd1, ncol=nb_sets , byrow = FALSE))
  colnames(sd2)<-sd_nm
  sd2<-data.frame(df[,1:6], sd2)
  return(sd2)}


path="C:/Users/mkaja/Dropbox/sHsp/mj_results/results/HDXMS/mia/all_data.csv"

names_states<- nm_states(path)
a5<-output_tp(path, states =names_states[c(5)], seq_match = F )
a1<-output_tp(path, states =names_states[c(1)], seq_match = F )


nb_ex1<-nb_exch_deut(a1)
av1<-ave_timepoint(a1)
sd1<-sd_timepoint(a1)
p1<-pv_timepoint(a1)
dav<-dif_ave(av1)
sd_ci<-sd_ci_timepoint(a1, alpha=0.05)
d1<-dim(av1)[2]


boxplot(c(a1[,7:9],a5[,7:9]))

#SD=sqrt(replicates)*(uper_limit-lower_limit)/t
# https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Confidence_Intervals/BS704_Confidence_Intervals5.html
# https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm

replicates=3
alpha=0.01

sd_ci<-sd_ci_timepoint(a1)

tvalue=abs(qt(alpha, replicates*2-2))
ts<-t.test(a1[15,7:9], a1[15,16:18], conf.level = c(1-alpha))
sd.ci<-(ts$conf.int[2]-ts$conf.int[1])/tvalue*sqrt(replicates)



#### anova tests

b1<-read.csv(path)
bt<-b1[which(b1$Deut.Time=="4.00s"),]

bt$nb_prot<-nb_exch_deut(bt)
bt$frac_deut<-as.numeric(bt$X..Deut)/bt$nb_prot



ggplot(bt) +
  aes(x = Experiment, y = frac_deut, color = Experiment) +
  geom_jitter() +
  theme(legend.position = "none")

boxplot(frac_deut ~ Experiment,
        data = bt
        )

kruskal.test(frac_deut ~ Experiment,
             data = bt
)
leveneTest(frac_deut ~ Protein.State,
           data = bt
)

res_aov <- aov(frac_deut ~ Experiment,
               data = bt
)

summary(res_aov)

par(mfrow = c(1, 2)) # combine plots

# 1. Homogeneity of variances
plot(res_aov, which = 3)

# 2. Normality
plot(res_aov, which = 2)

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(res_aov$residuals)

# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

shapiro.test(res_aov$residuals)



ggplot(bt) +
  aes(x = Experiment, y = frac_deut) +
  geom_boxplot()





#
plots_diff_tp(a1)

plot(av1[,7]/nb_ex1, pch=20)
points(av1[,10]/nb_ex1, col=5, pch=20)


plot(av1[,9]/nb_ex1, pch=20)
points(av1[,10]/nb_ex1, col=5, pch=20)


plot(av1[,8]/nb_ex1, col=2, pch=20)
points(av1[,9]/nb_ex1, col=4, pch=20)

plot(av1[,7]/nb_ex1, col=2, pch=20)
points(av1[,8]/nb_ex1, col=4, pch=20)


plot(av1[,7]/nb_ex1, av1[,8]/nb_ex1)


plot(av1[,7]/nb_ex1, av1[,10]/nb_ex1, pch=20)

plot(av1[,7]/nb_ex1, av1[,10]/nb_ex1, pch=20, xlab="", ylab="")
points(c(0,1),      c(0,1), lty=2, type="l")
points(c(0,1), c(0,1)+0.04, lty=2, type="l")
points(c(0,1), c(0,1)-0.04, lty=2, type="l")




head(av1) ###average differences (against first Protein State in the file)
da1<-dif_ave(av1)

CI_single(sd1[,7], 3)
CI_2pts(sd1[,7], sd1[,10],3)


nb_ex5<-nb_exch_deut(a5)


av5<-ave_timepoint(a5) sd5<-sd_timepoint(a5)

plot(av5[,7]/nb_ex5, av5[,8]/nb_ex5, pch=20, xlab="", ylab="") points(c(0,1),
                                                                      c(0,1), lty=2, type="l") points(c(0,1), c(0,1)+0.04, lty=2, type="l")
points(c(0,1), c(0,1)-0.04, lty=2, type="l")

plot(av5[,7]/nb_ex5, pch=20) points(av5[,8]/nb_ex5, col=2, pch=20)
points(av5[,9]/nb_ex5, col=4, pch=20) points(av5[,10]/nb_ex5, col=5, pch=20)

