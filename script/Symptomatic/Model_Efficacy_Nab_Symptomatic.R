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
library(ggrepel)

# Run from the `script` folder.
input_dir <- file.path("..", "processed data", "symptomatic efficacy")
output_dir <- file.path("Symptomatic", "Output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

### Import Tables For fitting Models
SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym<-read.csv(file.path(input_dir, "SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym.csv"))


### The Neutralisation ratio of vaccine to convalescence using reported neut titres
SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutConv)



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
LogisticEstimate=c("logk"=log(2.5),"C50"=log10(0.5))##InitialGuess

### Raw Efficacy Data Models - Fitting Logistic Model to Different combinations of Mean and SD
FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$MelbSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$MelbSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),LogisticEstimate),hessian=TRUE)


## CI raw efficacy logistic model data
FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,1)
FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,1)
FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,1)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,1)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,1)
FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1)

FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1]




####################################################################################################################################
######################Create Summary Table Of Fitted Estimate Titres
TableOfEstimatedTitresSym=data.frame("Method"=c(rep("Logistic_Raw",6)),
                                  "Mean"=c("CensoredEstimate","CensoredEstimate","CensoredEstimate","Reported","Reported","Reported"),
                                  "SD"=c(rep(c("byStudy_cens","MelbStudy_cens","Pooled_cens"),2)),
                                  "Estimate_EC50"=NA,
                                  "L_CI"=NA,
                                  "U_CI"=NA,
                                  "Estimate_k"=NA,
                                  "L_CI_k"=NA,
                                  "U_CI_k"=NA)



TableOfEstimatedTitresSym$Estimate_EC50[1]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,1)
TableOfEstimatedTitresSym$L_CI[1]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI[1]
TableOfEstimatedTitresSym$U_CI[1]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI[2]
TableOfEstimatedTitresSym$Estimate_k[1]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[1]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[1]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_CI_logk[2])

TableOfEstimatedTitresSym$Estimate_EC50[2]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,1)
TableOfEstimatedTitresSym$L_CI[2]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI[1]
TableOfEstimatedTitresSym$U_CI[2]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI[2]
TableOfEstimatedTitresSym$Estimate_k[2]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[2]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[2]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_CI_logk[2])

TableOfEstimatedTitresSym$Estimate_EC50[3]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,1)
TableOfEstimatedTitresSym$L_CI[3]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI[1]
TableOfEstimatedTitresSym$U_CI[3]=10^FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI[2]
TableOfEstimatedTitresSym$Estimate_k[3]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[3]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[3]=exp(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_CI_logk[2])

TableOfEstimatedTitresSym$Estimate_EC50[4]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,1)
TableOfEstimatedTitresSym$L_CI[4]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI[1]
TableOfEstimatedTitresSym$U_CI[4]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI[2]
TableOfEstimatedTitresSym$Estimate_k[4]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[4]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[4]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_CI_logk[2])

TableOfEstimatedTitresSym$Estimate_EC50[5]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,1)
TableOfEstimatedTitresSym$L_CI[5]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI[1]
TableOfEstimatedTitresSym$U_CI[5]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI[2]
TableOfEstimatedTitresSym$Estimate_k[5]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[5]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[5]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_CI_logk[2])

TableOfEstimatedTitresSym$Estimate_EC50[6]=10^tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1)
TableOfEstimatedTitresSym$L_CI[6]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI[1]
TableOfEstimatedTitresSym$U_CI[6]=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI[2]
TableOfEstimatedTitresSym$Estimate_k[6]=exp(tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1])
TableOfEstimatedTitresSym$L_CI_k[6]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI_logk[1])
TableOfEstimatedTitresSym$U_CI_k[6]=exp(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_CI_logk[2])


write.csv(TableOfEstimatedTitresSym,file=file.path(output_dir, "TableOfEstimatedTitresSym.csv"))



############################################################################################
####################Plotting Validation/Quality Of Fit Plots
###Logistic Model Raw Efficacy Data
FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SD,tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens$estimate,1))
FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$MelbSD,tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb$estimate,1))
FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_cens,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD,tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool$estimate,1))
FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SD,tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens$estimate,1))
FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$MelbSD,tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb$estimate,1))
FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD,tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1))


pdf(file.path(output_dir, "LogisticModel_RawEfficacy_QualityOfFit.pdf"),height=7,width=5)
par(mfrow=c(3,2))
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanCens_SDCens_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD by study" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanCens_SDMelb_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD Melb." )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanCens_SDPool_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD Pool" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanRept_SDCens_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD by study" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanRept_SDMelb_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD Melb." )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy,FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD Pool" )
abline(0,1,lty=2)
dev.off()





###### Manuscript Figures

NeutValue=seq(0.1,10,by=0.001)

SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$RatioReported_LB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported)-1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$RatioReported_UB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported)+1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$SEM)

Efficacy_Sym_Logistic<-NULL
Efficacy_Sym_Logistic_Raw<-NULL
for (i in 1:length(NeutValue)) {
  # Efficacy_Sym_Logistic[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[1],coef(FittedLogistic_MeanRept_SDMelb)[1],coef(FittedLogistic_MeanRept_SDMelb)[2])  
  Efficacy_Sym_Logistic_Raw[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),
                           SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[1],
                           tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                           tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}

LogisticModel_Sym_withPoolSD=data.frame("NeutRatio_Reported"=log10(NeutValue),"Efficacy"=Efficacy_Sym_Logistic_Raw)
LogisticModel_Sym_withPoolSD$Method<-rep("LogisticModel",length(NeutValue))

### Confidnece bounds:
Cov<-solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]

#################################################### Adding 95% Prediction Intervals #########################################################
grad1<-NULL
grad2<-NULL
Lower_Pred<-NULL
Upper_Pred<-NULL
G<-NULL

for (i in 1:length(NeutValue)) 
{
  f_temp <- function(p_temp) LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[1],p_temp[1],p_temp[2])
  grad1[i]<-numericGradient(f_temp, c(tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1)))[1]
  grad2[i]<-numericGradient(f_temp, c(tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1)))[2]
  G<-cbind(grad1[i],grad2[i])
  Lower_Pred[i]=Efficacy_Sym_Logistic_Raw[i]-1.96*sqrt(G%*%Cov%*%t(G))
  Upper_Pred[i]=Efficacy_Sym_Logistic_Raw[i]+1.96*sqrt(G%*%Cov%*%t(G))
}

LogisticModel_Sym_withPoolSD$Lower<-100*c(Lower_Pred)
LogisticModel_Sym_withPoolSD$Upper<-100*c(Upper_Pred) 
##########################################################################################################################################
Efficacy_fit_sym = NULL
for (i in 1:nrow(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym)){
  Efficacy_fit_sym[i] = LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported[i],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}
pearson_sym = cor.test(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy, Efficacy_fit_sym, method = "pearson")

rmse_sym = sqrt(mean(Efficacy_fit_sym*100 - SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Efficacy*100)^2)

ec50_sym_full <- 10^tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate, 1)
ec50_sym_full_label <- paste0("EC[50] == ", sprintf("%.2f", ec50_sym_full), "~'fold'")

Full_fit_sym<-ggplot(data=SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym, aes(y=100*Efficacy,x=(10^NeutRatio_Reported))) +
  #adding the bands
  geom_ribbon(data=LogisticModel_Sym_withPoolSD,aes(ymin=Lower, ymax=Upper), fill = "#D1E5F0", alpha = 0.55)+
  annotate("segment", x = 0.1, xend = ec50_sym_full, y = 50, yend = 50,
           linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  annotate("segment", x = ec50_sym_full, xend = ec50_sym_full, y = 0, yend = 50,
           linetype = "dashed", colour = "grey35", linewidth = 0.45) +
  annotate("text", x = ec50_sym_full * 1.18, y = 2.5, label = ec50_sym_full_label,
           hjust = 0, size = 4, parse = TRUE) +
  geom_point(shape=1) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper)) +
  geom_errorbarh(aes(xmin=RatioReported_LB,xmax=RatioReported_UB)) +
  scale_x_log10(lim=c(0.1,18),breaks=c(0.125,0.25,0.5,1,2,4,8),labels=c(0.125,0.25,0.5,1,2,4,8),
                expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(lim=c(0,110),breaks = seq(0, 100, by = 20),
                     expand = expansion(mult = c(0, 0.02))) +
  
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(),
        axis.title = element_text(size = 13, colour = "black"),
        axis.text = element_text(size = 12, colour = "black")) +
  geom_line(data=LogisticModel_Sym_withPoolSD, color="#2166AC", linewidth = 0.8) +
  ggrepel::geom_text_repel(aes(label=TechnicalName), size = 3, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.35, point.padding = 0.2, force = 5, seed = 1234) +
  xlab("Normalised neutralisation level (fold of convalescent)") +
  ylab("Protective efficacy against symptomatic infection (%)")

pdf(file.path(output_dir, "Full_fit_sym.pdf"),height=5,width=5.5)
print(Full_fit_sym)
dev.off()




####This creates EC50 estimates for each cohort 

fixlogk_Sym=tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model<-nlm(function(p){FittingLogistic_Raw(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)]),fixlogk_Sym,p[(1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)))],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],
                                                                                          SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)])},rep(LogisticEstimate[2],sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))),hessian=TRUE)

FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model_CI<- matrix(c(-tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)))*1.96+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont))),tail(sqrt(diag(solve(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)))*1.96+tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$NumCont)))),ncol=2)

FittingEachStudySep_LogisticRaw_ReptMean_PooledSD=data.frame("Study"=SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Study,"EC50"=10^tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,length(SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$Study)),"CI_L"=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,1],"CI_U"=10^FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,2])

write.csv(FittingEachStudySep_LogisticRaw_ReptMean_PooledSD,file=file.path(output_dir, "FittingEachStudySep_LogisticRaw_ReptMean_PooledSD.csv"))

#####Protective Titres by study
SuppTable_TableOfSDperStudy<-read.csv(file.path(input_dir, "SuppTable_TableOfSDperStudy.csv"))
SuppTable_TableOfSDperStudy$ProtectiveTitre<-SuppTable_TableOfSDperStudy$NeutConv*TableOfEstimatedTitresSym$Estimate_EC50[TableOfEstimatedTitresSym$Method=="Logistic_Raw" & TableOfEstimatedTitresSym$Mean=="Reported" & TableOfEstimatedTitresSym$SD=="Pooled_cens"]
SuppTable_TableOfSDperStudy$ProtectiveTitre<-round(SuppTable_TableOfSDperStudy$ProtectiveTitre,0)
write.csv(SuppTable_TableOfSDperStudy,file.path(output_dir, "SuppTable_TableOfSDperStudy_withProtectiveTitres.csv"))




#### Predict efficacy of fractional doses

points_Sym<- read.csv(file.path(input_dir, "PooledData_Doses.csv"))
points_Sym['type'] <- 'Alternative doses'
points_Sym[points_Sym$dose==1,'type'] <- 'Status quo dose'

platform_palette <- c(
  "mRNA" = "#8C510A",
  "Protein subunit" = "#01665E",
  "Viral vector" = "#5E3C99",
  "Inactivated" = "#D73027",
  "DNA" = "#2C7BB6",
  "Convalescent" = "#4D4D4D"
)

Efficacy_pred_Sym <- NULL
for (i in 1:length(points_Sym$neutralization)){
  Efficacy_pred_Sym[i]=LogisticModel_PercentUninfected(points_Sym$neutralization[i],
                                                       SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym$PooledSD[1],
                                                       tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,2)[1],
                                                       tail(FittedLogistic_Sym_RawEfficacy_MeanRept_SDPool$estimate,1))  
}

points_Sym$efficacy <- Efficacy_pred_Sym

Frac_pred_sym <- ggplot() +
  geom_ribbon(data = LogisticModel_Sym_withPoolSD, aes(ymin = Lower, ymax = Upper, x = 10^NeutRatio_Reported), fill = "grey 90") +
  geom_line(data = LogisticModel_Sym_withPoolSD, aes(x=10^NeutRatio_Reported,y=100*Efficacy), colour = "grey20", linewidth=0.8) +
  geom_point(data=points_Sym,
             aes(x=10^neutralization,y=100*efficacy,color=vaccine_type,shape=type,size=dose)) +
  scale_color_manual(values = platform_palette, name = "Vaccine platform") +
  scale_shape_manual(name="",values=c(16, 10)) +
  scale_size_continuous(range=c(3.5,6)) +
  ylab("Protective efficacy against symptomatic infection (%)") + xlab("Normalised neutralisation level (fold of convalescent)") +
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

pdf(file.path(output_dir, "Frac_pred_sym.pdf"),height=5.8,width=5.8)
print(Frac_pred_sym)
dev.off()
