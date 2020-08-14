library(survminer)
set.seed(258)

#change the order of column/patients of clinical data following the variable 'info'
cli = c_clinical; rownames(cli) = c_clinical$PATIENT_ID
grp=sub_grp; grp=as.data.frame(grp)
cli_surv =cli_surv[rownames(grp),]
survData<-cbind(cli_surv$OS_MONTHS, cli_surv$OS_STATUS=="DECEASED")

#create survData
rownames(survData)<- rownames(cli_surv)
colnames(survData)<-c("time", "status")
survData <- as.data.frame(survData)

coxFit1 <- coxph(
  Surv(time, status) ~ as.factor(sub_grp),
  data = survData,
  ties = "exact"
)
pcox=round(summary(coxFit1)$logtest[3],8);pcox #p-value = 1.498e-05  

mfit <- survfit(Surv(time, status == 1) ~ as.factor(sub_grp), data = survData)

#visualization of survival rate between two groups
ggsurvplot(mfit, size=1,
           linetype = "strata",
           risk.table = FALSE, fun ="pct", risk.table.col = "strata", break.x.by = 50,
           xlab = "Time in months",
           legend = "bottom",
           legend.title = "Group", legend.labs = c("Group 1","Group 2"),
           conf.int = TRUE, pval = paste("P-value",pcox,sep=" = " ), xlim = c(0,200),
           palette = c("black", "green"))