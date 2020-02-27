library(dplyr)
library(ggplot2)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library("ggpubr")
library(ggsci)

#prepare data
a =as.data.frame(sub_grp)
colnames(a) <- c('Groups')
a$lymph=c_clinical$LYMPH_NODES_EXAMINED_POSITIVE
a$npi=c_clinical$NPI
a$stage=c_clinical$stage
head(a)
#visualize
set.seed(216)
#lymph nodes
my_comparisons <- list(c("1", "2"))
p=ggboxplot(a, x = "Groups", y = "lymph", 
          color = "Groups", palette = "npg",
          ylab = "Number of positive lymph nodes", xlab = "Groups") + border("black")
ggpar(p,legend="right",legend.title = "group")
p + stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(aes(label=paste0("Wilcoxon, P-value = ",..p.format..)),label.x = 1.3,label.y = 52)

#NPI
p1<-ggboxplot(a, x = "Groups", y = "npi", 
                     color = "Groups", palette = "npg",
                     ylab = "Nottingham Prognostic Index", xlab = "Groups",
                    bxp.errorbar=TRUE) + border("black")
ggpar(p1,legend="right", legend.title = "group")
p1 + stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(aes(label=paste0("Wilcoxon, P-value = ",..p.format..)),label.x = 1.3,label.y = 8)

#stage
a_st = a %>% group_by(Groups, stage) %>% summarise(patient.num = n()) ## Create a column "patient.num" which is number of patients in each combination of group and stage
a_st$Groups=as.character(a_st$Groups)
a_st$stage=as.character(a_st$stage)
a_st=na.omit(a_st)
p2=ggplot(a_st, aes(x = Groups, y = patient.num, fill = stage)) +   #stacked barplot
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format())
p2 + ylab("Percent of patients")

library(MASS)       # load the MASS package 
#test chi-square
set.seed(357)
achitest = data.frame(a$Groups, a$stage)
achitest=na.omit(achitest)
colnames(achitest) = c("Groups", "stage")
achitest$Groups = as.character(achitest$Groups)
achitest$stage = as.character(achitest$stage)
tbl = table(achitest$Groups, achitest$stage)
chisq.test(tbl)