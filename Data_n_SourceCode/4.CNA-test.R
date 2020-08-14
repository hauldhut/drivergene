#library
devtools::install_github("kassambara/ggpubr")
library("ggplot2")
library("ggpubr")
library("ggsci")

#prepare data
a = info
cli_feature = c_clinical; rownames(cli_feature) = cli_feature$PATIENT_ID; cli_feature = cli_feature[rownames(info),]
a$lymph=cli_feature$LYMPH_NODES_EXAMINED_POSITIVE
a$npi=cli_feature$NPI
a$stage=cli_feature$stage
head(a)

#comparison
install.packages("table1"); library(table1)
install.packages("compareGroups"); library("compareGroups")
#define specifically type of data
a$lymph = as.numeric(a$lymph); a$npi = as.numeric(a$npi)
a$stage = as.character(a$stage) 
#start to perform comparisons
des=compareGroups::createTable(compareGroups::compareGroups(groups ~ ., data = a, method = NA))
#save the results as xls file
compareGroups::export2xls(des, file = "tableSTAT.xlsx", header.labels = c(p.overall = "p-value"))

#visualize
set.seed(216)
#lymph nodes
p=ggboxplot(a, x = "groups", y = "lymph", 
          color = "groups", palette = c("black", "green"),
          ylab = "Number of positive lymph nodes", xlab = "Groups",
          title = "Wilcoxon, P-value = 0.031") + border("black")
ggpar(p,legend="right",legend.title = "Groups")

#NPI
p1=ggboxplot(a, x = "groups", y = "npi", 
            color = "groups", palette = c("black", "green"),
            ylab = "Nottingham prognostic index", xlab = "Groups",
            title = "Wilcoxon, P-value < 0.001") + border("black")
ggpar(p1,legend="right",legend.title = "Groups")

#stage
a_st = a %>% group_by(groups, stage) %>% summarise(patient.num = n()) ## Create a column "patient.num" which is number of patients in each combination of group and stage
a_st$Groups=as.character(a_st$groups)
a_st$stage=as.character(a_st$stage)
a_st=na.omit(a_st)
p2=ggplot(a_st, aes(x = Groups, y = patient.num, fill = stage)) +   #stacked barplot
  geom_bar(position = "fill", stat = "identity") + theme(panel.background = element_blank()) +
  ggtitle("Chisq, P-value = 0.016") +
  scale_y_continuous(labels = scales::percent_format())
p2 + xlab("Groups") + ylab("Percent of patients")
