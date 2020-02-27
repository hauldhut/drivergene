#library
if(!require(survival)) install.packages('survival')
if(!require(Hmisc)) install.packages('Hmisc')
library('survival')
library('Hmisc')
library(survminer) 

set.seed(258)
survData<-cbind(c_clinical$OS_MONTHS, c_clinical$OS_STATUS=="DECEASED")

rownames(survData)<-rownames(c_clinical)
colnames(survData)<-c("time", "status")
survData <- as.data.frame(survData)

coxFit1 <- coxph(
  Surv(time, status) ~ as.factor(sub_grp),
  data = survData,
  ties = "exact"
)
pcox=round(summary(coxFit1)$coefficients[5],5) #p-value

mfit <- survfit(Surv(time, status == 1) ~ as.factor(sub_grp), data = survData)

ggsurvplot(mfit, size=1,
           linetype = "strata",
           risk.table = TRUE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in months",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("Group 1","Group 2"),
           conf.int = TRUE, pval = paste("P-value",pcox,sep=" = " ), ggtheme = theme_gray(),xlim = c(0,200)) 