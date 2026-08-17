library(ggplot2)
library(ggpubr)
library(ggsci)
library(scales)
library(lemon)
library(shades)
library(plyr)
library(car)
library(maxLik) #To compute numerical gradient as a function of the parameters 
library(tidyverse)
library(meta)
library(ggrepel)

# Run from the `script` folder.
input_dir <- file.path("..", "processed data", "asymptomatic efficacy")
output_dir <- file.path("Asymptomatic", "Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

### Import Tables For fitting Models
Summary_Efficacy_NeutRatio_SD_SEM_Asym<-read.csv(file.path(input_dir, "Summary_Efficacy_NeutRatio_SD_SEM_Asym.csv")) ##For threshold and logistic


### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported=log10(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutMean/Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutConv)

Summary_Efficacy_NeutRatio_SD_SEM_Asym$SeSD[3]=mean((Summary_Efficacy_NeutRatio_SD_SEM_Asym %>% dplyr::filter(!is.na(SeSD)))$SeSD)

Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool=metagen(TE = SD, seTE = SeSD, studlab = Study, 
        data = Summary_Efficacy_NeutRatio_SD_SEM_Asym,
        n.e = NumberIndividuals_Vaccine, level.ci = 0.95, hakn = TRUE)$TE.random
##########################################################################################
##################### The Models

####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  Output<-NULL
  
  if (length(C50)==1) {
    C50=rep(C50,length(mu_titre))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(mu_titre))
  }
  
  for (i in 1:length(mu_titre)) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}


### Logistic model for Raw Efficacy Counts
FittingLogistic_Raw<-function(logRisk0,logk,C50,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector){
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(N_C))
  }

  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  if (length(C50)==1) {
    C50=rep(C50,length(N_C))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(N_C))
  }
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk[i],C50[i]))))
  }
  LL
}



##########################################################################################
##Fitting Models
LogisticEstimate_Asym = c("logk"=log(1.5),"C50"=log10(2.5))##InitialGuess

### Raw Efficacy Data Models - Fitting Logistic Model to Different SD (byStudy SD & pooled SD)


FittedLogistic_Asym_RawEfficacy_MeanRept_SD<-nlm(function(p)
  {FittingLogistic_Raw(p[1:sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))],
                       p[sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))+1],
                       p[sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))+2],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                       Summary_Efficacy_NeutRatio_SD_SEM_Asym$SD[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)])},
  c(log(Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]/
          Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]),
    LogisticEstimate_Asym),hessian=TRUE)

FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool<-nlm(function(p)
{FittingLogistic_Raw(p[1:sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))],
                     p[sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))+1],
                     p[sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))+2],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
                     Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)])},
c(log(Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]/
        Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]),
  LogisticEstimate_Asym),hessian=TRUE)

## CI raw efficacy logistic model data
FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI <- tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$hessian))),1)*
  1.96*c(-1,1)+tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,1)
FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI_logk <- tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$hessian))),2)[1]*
  1.96*c(-1,1)+tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,2)[1]

FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI<- tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$hessian))),1)*
  1.96*c(-1,1)+tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1)
FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI_logk <- tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$hessian))),2)[1]*
  1.96*c(-1,1)+tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1]

####################################################################################################################################
######################Create Summary Table Of Fitted Estimate Titres
TableOfEstimatedTitresAsym=data.frame("Method"=c("Logistic_Raw","Logistic_Raw"),
                                  "Mean"=c("Reported","Reported"),
                                  "SD"=c("byStudy_cens","Pooled_cens"),
                                  "Estimate_EC50"=NA,
                                  "L_CI"=NA,
                                  "U_CI"=NA,
                                  "Estimate_k"=NA,
                                  "L_CI_k"=NA,
                                  "U_CI_k"=NA)


TableOfEstimatedTitresAsym$Estimate_EC50[1]=10^tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,1)
TableOfEstimatedTitresAsym$L_CI[1]=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI[1]
TableOfEstimatedTitresAsym$U_CI[1]=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI[2]
TableOfEstimatedTitresAsym$Estimate_k[1]=exp(tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,2)[1])
TableOfEstimatedTitresAsym$L_CI_k[1]=exp(FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI_logk[1])
TableOfEstimatedTitresAsym$U_CI_k[1]=exp(FittedLogistic_Asym_RawEfficacy_MeanRept_SD_CI_logk[2])

TableOfEstimatedTitresAsym$Estimate_EC50[2]=10^tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1)
TableOfEstimatedTitresAsym$L_CI[2]=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI[1]
TableOfEstimatedTitresAsym$U_CI[2]=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI[2]
TableOfEstimatedTitresAsym$Estimate_k[2]=exp(tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1])
TableOfEstimatedTitresAsym$L_CI_k[2]=exp(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI_logk[1])
TableOfEstimatedTitresAsym$U_CI_k[2]=exp(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_CI_logk[2])

write.csv(TableOfEstimatedTitresAsym,file=file.path(output_dir, "TableOfEstimatedTitresAsym.csv"))



############################################################################################
####################Plotting Validation/Quality Of Fit Plots
###Logistic Model Raw Efficacy Data

FittedLogistic_Asym_RawEfficacy_MeanRept_SD_ModelOutput<-
  LogisticModel_PercentUninfected(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported,
                                  Summary_Efficacy_NeutRatio_SD_SEM_Asym$SD,
                                  tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,2)[1],
                                  tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SD$estimate,1))

FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_ModelOutput<-
  LogisticModel_PercentUninfected(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported,
                                  Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool,
                                  tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                  tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1))

pdf(file.path(output_dir, "LogisticModel_Asym_RawEfficacy_QualityOfFit.pdf"))
par(mfrow=c(1,2), pty = 's')
plot(Summary_Efficacy_NeutRatio_SD_SEM_Asym$Efficacy,FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_ModelOutput,
     xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD by study" )
abline(0,1,lty=2)
plot(Summary_Efficacy_NeutRatio_SD_SEM_Asym$Efficacy,FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_ModelOutput,
     xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, Pooled SD " )
abline(0,1,lty=2)
dev.off()


###### Manuscript Figures

NeutValue=seq(0.1,10,by=0.001)

Summary_Efficacy_NeutRatio_SD_SEM_Asym$RatioReported_LB=
  10^((Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported)
      -1.96*Summary_Efficacy_NeutRatio_SD_SEM_Asym$SEM)

Summary_Efficacy_NeutRatio_SD_SEM_Asym$RatioReported_UB=
  10^((Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported)
     +1.96*Summary_Efficacy_NeutRatio_SD_SEM_Asym$SEM)

Efficacy_Asym_Logistic<-NULL
Efficacy_Asym_Logistic_Raw<-NULL
for (i in 1:length(NeutValue)) {
    Efficacy_Asym_Logistic_Raw[i]=
    LogisticModel_PercentUninfected(log10(NeutValue[i]),
    Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool[1],
    tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
    tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}

LogisticModel_Asym_withPooledSD = data.frame("NeutRatio_Reported"=log10(NeutValue),"Efficacy"=Efficacy_Asym_Logistic_Raw)
LogisticModel_Asym_withPooledSD$Method<-rep("LogisticModel",length(NeutValue))


### Confidnece bounds:
Cov_asym<-solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$hessian)[7:8,7:8]


#################################################### Adding 95% Prediction Intervals #########################################################
grad1_Asym<-NULL
grad2_Asym<-NULL
Lower_Pred_Asym<-NULL
Upper_Pred_Asym<-NULL
G_asym<-NULL

for (i in 1:length(NeutValue)) 
{
  f_temp <- function(p_temp) 
    LogisticModel_PercentUninfected(log10(NeutValue[i]),Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool[1],p_temp[1],p_temp[2])
  grad1_Asym[i]<-numericGradient(f_temp, c(tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1)))[1]
  grad2_Asym[i]<-numericGradient(f_temp, c(tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1)))[2]
  G_asym<-cbind(grad1_Asym[i],grad2_Asym[i])
  Lower_Pred_Asym[i]=Efficacy_Asym_Logistic_Raw[i]-1.96*sqrt(G_asym%*%Cov_asym%*%t(G_asym))
  Upper_Pred_Asym[i]=Efficacy_Asym_Logistic_Raw[i]+1.96*sqrt(G_asym%*%Cov_asym%*%t(G_asym))
}

LogisticModel_Asym_withPooledSD$Lower<-100*c(Lower_Pred_Asym)
LogisticModel_Asym_withPooledSD$Upper<-100*c(Upper_Pred_Asym) 
##########################################################################################################################################
Efficacy_fit_asym = NULL
for (i in 1:nrow(Summary_Efficacy_NeutRatio_SD_SEM_Asym)){
  Efficacy_fit_asym[i] = LogisticModel_PercentUninfected(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported[i],
                                                        Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool[1],
                                                        tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                                        tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}
pearson_asym = cor.test(Summary_Efficacy_NeutRatio_SD_SEM_Asym$Efficacy, Efficacy_fit_asym, method = "pearson")

rmse_asym = sqrt(mean(Efficacy_fit_asym*100 - Summary_Efficacy_NeutRatio_SD_SEM_Asym$Efficacy*100)^2)

ec50_asym_full <- 10^tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate, 1)
ec50_asym_full_label <- paste0("EC[50] == ", sprintf("%.2f", ec50_asym_full), "~'fold'")

Full_fit_asym<-ggplot(data=Summary_Efficacy_NeutRatio_SD_SEM_Asym, aes(y=100*Efficacy,x=(10^NeutRatio_Reported))) +
  #adding the bands
  geom_ribbon(data=LogisticModel_Asym_withPooledSD,aes(ymin=Lower, ymax=Upper), fill = "#D1E5F0", alpha = 0.55)+
  annotate("segment", x = 0.1, xend = ec50_asym_full, y = 50, yend = 50,
           linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  annotate("segment", x = ec50_asym_full, xend = ec50_asym_full, y = 0, yend = 50,
           linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  annotate("text", x = ec50_asym_full * 1.08, y = 2.5, label = ec50_asym_full_label,
           hjust = 0, size = 4, parse = TRUE) +
  geom_point(shape=1) +
  geom_errorbar(aes(ymin=ifelse(Lower<0,0,Lower),ymax=Upper)) +
  geom_errorbarh(aes(xmin=RatioReported_LB,xmax=RatioReported_UB)) +
  scale_x_log10(lim=c(0.1,13),breaks=c(0.125,0.25,0.5,1,2,4,8),labels=c(0.125,0.25,0.5,1,2,4,8),
                expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20),
                     expand = expansion(mult = c(0, 0.02))) +
  
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        axis.title = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 12, colour = "black")) +
  geom_line(data=LogisticModel_Asym_withPooledSD,color="#2166AC", linewidth = 0.8) +
  ggrepel::geom_text_repel(aes(label=Vaccine), size = 3, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.2, point.padding = 0.15) +
  xlab("Normalised neutralisation level (fold of convalescent)") +
  ylab("Protective efficacy against asymptomatic infection (%)")

pdf(file.path(output_dir, "Full_fit_asym.pdf"),height=5,width=5.5)
print(Full_fit_asym)
dev.off()




####This creates EC50 estimates for each cohort 

fixlogk_Asym = tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model<-nlm(function(p){FittingLogistic_Raw(
  log(Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]/
        Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)]),
  fixlogk_Asym,p[(1:sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)))],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfCont[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$InfVac[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$NeutRatio_Reported[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)],
  Summary_Efficacy_NeutRatio_SD_SEM_Asym$SD[!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)])},
  rep(LogisticEstimate_Asym[2],sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))),hessian=TRUE)

FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model_CI<- matrix(
  c(-tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),
          sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)))*1.96+
      tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,
           sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont))),
    tail(sqrt(diag(solve(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),
         sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)))*1.96+
      tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,
           sum(!is.na(Summary_Efficacy_NeutRatio_SD_SEM_Asym$NumCont)))),ncol=2)

FittingEachStudySep_Asym_LogisticRaw_ReptMean_SDPool=data.frame(
  "Study"=Summary_Efficacy_NeutRatio_SD_SEM_Asym$Study,
  "EC50"=10^tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,
                 length(Summary_Efficacy_NeutRatio_SD_SEM_Asym$Study)),
  "CI_L"=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,1],
  "CI_U"=10^FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,2])

write.csv(FittingEachStudySep_Asym_LogisticRaw_ReptMean_SDPool,file=file.path(output_dir, "FittingEachStudySep_Asym_LogisticRaw_ReptMean_SDPool.csv"))

##### 50% Protective Titres by study
SuppTable_Asym_TableOfSDperStudy<-read.csv(file.path(input_dir, "SuppTable_Asym_TableOfSDperStudy.csv"))
SuppTable_Asym_TableOfSDperStudy$ProtectiveTitre<-SuppTable_Asym_TableOfSDperStudy$NeutConv*TableOfEstimatedTitresAsym$Estimate_EC50[TableOfEstimatedTitresAsym$Method=="Logistic_Raw" & TableOfEstimatedTitresAsym$Mean=="Reported" & TableOfEstimatedTitresAsym$SD=="Pooled_cens"]
SuppTable_Asym_TableOfSDperStudy$ProtectiveTitre<-round(SuppTable_Asym_TableOfSDperStudy$ProtectiveTitre,0)
write.csv(SuppTable_Asym_TableOfSDperStudy,file.path(output_dir, "SuppTable_Asym_TableOfSDperStudy_withProtectiveTitres.csv"))



#### Predict efficacy of fractional doses

points_asym<- read.csv(file.path(input_dir, "PooledData_Doses.csv"))
points_asym['type'] <- 'Alternative doses'
points_asym[points_asym$dose==1,'type'] <- 'Status quo dose'

platform_palette <- c(
  "mRNA" = "#8C510A",
  "Protein subunit" = "#01665E",
  "Viral vector" = "#5E3C99",
  "Inactivated" = "#D73027",
  "DNA" = "#2C7BB6",
  "Convalescent" = "#4D4D4D"
)

Efficacy_pred_asym <- NULL
for (i in 1:length(points_asym$neutralization)){
  Efficacy_pred_asym[i]=LogisticModel_PercentUninfected(points_asym$neutralization[i],
                                              Summary_Efficacy_NeutRatio_SD_SEM_Asym$SDPool[1],
                                              tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                              tail(FittedLogistic_Asym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}

points_asym$efficacy <- Efficacy_pred_asym

Frac_pred_asym <- ggplot() +
  geom_ribbon(data = LogisticModel_Asym_withPooledSD, aes(ymin = Lower, ymax = Upper, x = 10^NeutRatio_Reported), fill = "grey 90") +
  geom_line(data = LogisticModel_Asym_withPooledSD, aes(x=10^NeutRatio_Reported,y=100*Efficacy), colour = "grey20", linewidth=0.8) +
  geom_point(data=points_asym,
             aes(x=10^neutralization,y=100*efficacy,color=vaccine_type,shape=type,size=dose)) +
  scale_color_manual(values = platform_palette, name = "Vaccine platform") +
  scale_shape_manual(name="",values=c(16, 10)) +
  scale_size_continuous(range=c(3.5,6)) +
  ylab("Protective efficacy against asymptomatic infection (%)") + xlab("Normalised neutralisation level (fold of convalescent)") +
  theme(panel.grid.major=element_line(colour= "grey 92", linewidth = 0.28),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        aspect.ratio = 1,
        legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        legend.text = element_text(size=12),
        legend.key.size = grid::unit(0.65, "cm"),
        axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        plot.margin = margin(8, 12, 8, 12)) +
  guides(color = "none", size = "none", shape=guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0, 0.025))) +
  scale_x_log10(lim=c(0.1,10),breaks=c(0.1, 1, 10),labels=c(0.1, 1.0, 10.0), expand = c(0, 0))

pdf(file.path(output_dir, "Frac_pred_asym.pdf"),height=5.8,width=5.8)
print(Frac_pred_asym)
dev.off()
