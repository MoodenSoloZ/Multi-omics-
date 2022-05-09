survive_result<-function(clinicali,label){
  clinicali <- cbind(clinicali,data.frame(cluster=label))
  clinicali$Death <- as.numeric(as.character(clinicali$Death))
  coxFit <- coxph(Surv(time = Time, event = Death) ~ as.factor(cluster), data = clinicali, ties = "exact")
  coxresults<-summary(coxFit)
  fit1<-survfit(Surv(Time, Death) ~ as.factor(cluster), data = clinicali)
  g<-ggsurvplot(fit1,
                #survfit(coxFit, newdata = data.frame(cluster=c(1,2,3,4))),
                data =clinicali,
                palette = seq(1,length(summary(as.factor(label)))),
                conf.int =FALSE,
                xlab ="Time(Days)", 
                legend.title = "Survival curve",
                legend.labs=seq(1,length(summary(as.factor(label)))) 
  )
  g$plot <- g$plot +annotate("text", x = 5000, y = 1.0, 
                             label = paste("Pvalue=",signif(coxresults$logtest[3],3)),
                             cex=5, vjust=0, hjust = 1.1, fontface=2)
  print(coxresults$logtest[3])
  g
}


