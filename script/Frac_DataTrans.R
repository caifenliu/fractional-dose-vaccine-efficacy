library(readxl)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(scales)
library(lemon)
library(ggsci)
library(meta)

# Run from the `script` folder.
data_dir <- file.path("..", "processed data", "dose neutralisation")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

data_raw = read_csv(file.path(data_dir, "data_raw.csv"))
data_raw = data_raw%>%drop_na(nab, conv)
data_raw$logTiter_vaccine = log10(data_raw$nab)
data_raw$logTiter_conv = log10(data_raw$conv)
data_raw$sd_vaccine = sqrt(data_raw$N)*(log10(data_raw$nab_upper)-log10(data_raw$nab_lower))/3.92
data_raw$sd_conv = sqrt(data_raw$conv_N)*(log10(data_raw$conv_upper)-log10(data_raw$conv_lower))/3.92
data_raw$se_vaccine = (log10(data_raw$nab_upper)-log10(data_raw$nab_lower))/3.92
data_raw$se_conv = (log10(data_raw$conv_upper)-log10(data_raw$conv_lower))/3.92
data_raw$log_ratio_to_conv = log10(data_raw$nab)-log10(data_raw$conv)
data_raw$ratio_to_conv = 10^data_raw$log_ratio_to_conv
data_raw$se_sd_conv = data_raw$sd_conv / sqrt(2 * (data_raw$conv_N - 1))
data_raw$se_log_ratio_to_conv = sqrt(data_raw$se_vaccine^2 + data_raw$se_conv^2)

bootstrap_response_sd <- function(mu, response_sd, N, n_boot = 2000) {
  if (is.na(mu) || is.na(response_sd) || is.na(N) || N <= 1 || response_sd <= 0) {
    return(c(mean_sd_vaccine = NA, se_sd_vaccine = NA))
  }

  log_ratios <- matrix(
    rnorm(N * n_boot, mean = mu, sd = response_sd),
    nrow = N,
    ncol = n_boot
  )
  hat_sd_log_ratios <- apply(log_ratios, 2, sd)

  c(
    mean_sd_vaccine = mean(hat_sd_log_ratios),
    se_sd_vaccine = sd(hat_sd_log_ratios)
  )
}
  

transNewVars = function(data){
  
  set.seed(1234)
  used = data %>% arrange(study_id, measured_day)
  
  
  ## Bootstrap on mean neutralising titre (log ratios)
  
  temp = t(sapply(1:nrow(used), function(j){
    
    mu = used$log_ratio_to_conv[j]
    sd = used$sd_vaccine[j]
    N = used$N[j]
    
    bootstrap_response_sd(mu, sd, N)
    
  }))
  
  
  colnames(temp) = c("mean_sd_vaccine", "se_sd_vaccine")
  used = cbind(used, temp)
  used
  
  pool_sd_vaccine <- metagen(TE = sd_vaccine, seTE = se_sd_vaccine, studlab = study, 
                             data = used %>% dplyr::filter(!is.na(sd_vaccine) & sd_vaccine != 0 & !is.na(se_sd_vaccine) & se_sd_vaccine != 0) ,
                             n.e = N, level.ci = 0.95, hakn = TRUE)$TE.random

  pool_sd_conv <- metagen(TE = sd_conv, seTE = se_sd_conv, studlab = study,
                          data = used %>% dplyr::filter(!is.na(sd_conv) & sd_conv != 0 & !is.na(se_sd_conv) & se_sd_conv != 0),
                          n.e = conv_N, level.ci = 0.95, hakn = TRUE)$TE.random
  
  used$sd_vaccine = ifelse(is.na(used$sd_vaccine)|used$sd_vaccine == 0, pool_sd_vaccine, used$sd_vaccine)
  used$sd_conv = ifelse(is.na(used$sd_conv)|used$sd_conv == 0, pool_sd_conv, used$sd_conv)
  used$se_vaccine = used$sd_vaccine/sqrt(used$N)
  used$se_conv = ifelse(!is.na(used$conv_N) & used$conv_N > 0, used$sd_conv/sqrt(used$conv_N), used$se_conv)
  used$se_log_ratio_to_conv = sqrt(used$se_vaccine^2 + used$se_conv^2)
  used$log_ratio_lower = used$log_ratio_to_conv - 1.96*used$se_log_ratio_to_conv
  used$log_ratio_upper = used$log_ratio_to_conv + 1.96*used$se_log_ratio_to_conv
  
  
  set.seed(1234)
  for(j in 1:nrow(used)){
    
    if (is.na(used$mean_sd_vaccine[j]) | used$mean_sd_vaccine[j] == 0) {
      
      mu = used$log_ratio_to_conv[j]
      sd = used$sd_vaccine[j]
      N = used$N[j]
      
      sd_summary = bootstrap_response_sd(mu, sd, N)
      
      used$mean_sd_vaccine[j] = sd_summary["mean_sd_vaccine"]
      used$se_sd_vaccine[j] = sd_summary["se_sd_vaccine"]
      
    }} 
  
  
  test = do.call('cbind', used)
  test = as.data.frame(used)
  
  return(test)
  
}
   


data_trans = transNewVars(data = data_raw)
write.csv(data_trans, file = file.path(data_dir, "data_trans.csv"))




# Get pooled estimates by vaccine type and fractional dose   

pool_log_ratio_frac_day = function(data, day){
  output = metagen(TE = log_ratio_to_conv, seTE = se_log_ratio_to_conv, studlab = study, 
          data = data %>% dplyr::filter(measured_day == day),
          n.e = N, level.ci = 0.95, hakn = TRUE)
  return(output)
}


## mRNA, fraction = 0.25, 0.33, 0.4, 0.5, 0.67, 1

data_mRNA = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="mRNA")

write.csv(data_mRNA, file = file.path(data_dir, "data_mRNA.csv"))



mRNA_frac_less_0.5 = data_mRNA %>%
  dplyr::filter(dose_frac < 0.5)

mRNA_frac_greater_0.5 = data_mRNA %>%
  dplyr::filter(dose_frac >= 0.5 & dose_frac < 1)

mRNA_full = data_mRNA %>%
  dplyr::filter(dose_frac == 1)



pool_log_ratio_mRNA_frac_less_0.5_day0 = pool_log_ratio_frac_day(mRNA_frac_less_0.5, 0)
pool_log_ratio_mRNA_frac_less_0.5_day14 = pool_log_ratio_frac_day(mRNA_frac_less_0.5, 14)
pool_log_ratio_mRNA_frac_less_0.5_day21 = pool_log_ratio_frac_day(mRNA_frac_less_0.5, 21)
pool_log_ratio_mRNA_frac_less_0.5_day28 = pool_log_ratio_frac_day(mRNA_frac_less_0.5, 28)
pool_log_ratio_mRNA_frac_less_0.5_day35 = pool_log_ratio_frac_day(mRNA_frac_less_0.5, 35)
pool_log_ratio_mRNA_frac_less_0.5_day42 <- pool_log_ratio_frac_day(mRNA_frac_less_0.5, 42)
pool_log_ratio_mRNA_frac_less_0.5_day56 <- pool_log_ratio_frac_day(mRNA_frac_less_0.5, 56)


pool_log_ratio_mRNA_frac_greater_0.5_day0 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 0)
pool_log_ratio_mRNA_frac_greater_0.5_day21 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 21)
pool_log_ratio_mRNA_frac_greater_0.5_day28 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 28)
pool_log_ratio_mRNA_frac_greater_0.5_day35 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 35)
pool_log_ratio_mRNA_frac_greater_0.5_day42 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 42)
pool_log_ratio_mRNA_frac_greater_0.5_day56 = pool_log_ratio_frac_day(mRNA_frac_greater_0.5, 56)


pool_log_ratio_mRNA_full_day0 = pool_log_ratio_frac_day(mRNA_full, 0)
pool_log_ratio_mRNA_full_day14 = pool_log_ratio_frac_day(mRNA_full, 14)
pool_log_ratio_mRNA_full_day21 = pool_log_ratio_frac_day(mRNA_full, 21)
pool_log_ratio_mRNA_full_day28 = pool_log_ratio_frac_day(mRNA_full, 28)
pool_log_ratio_mRNA_full_day35 = pool_log_ratio_frac_day(mRNA_full, 35)
pool_log_ratio_mRNA_full_day42 = pool_log_ratio_frac_day(mRNA_full, 42)
pool_log_ratio_mRNA_full_day56 = pool_log_ratio_frac_day(mRNA_full, 56)



## Inactivated vaccine, fraction =  0.5, 1

data_inact = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="Inactivated vaccine")

write.csv(data_inact, file = file.path(data_dir, "data_inact.csv"))


inact_frac_less_0.5 = data_inact %>%
  dplyr::filter(dose_frac < 0.5)

inact_frac_greater_0.5 = data_inact %>%
  dplyr::filter(dose_frac >= 0.5 & dose_frac < 1)

inact_full = data_inact %>%
  dplyr::filter(dose_frac == 1)

inact_frac_2 = data_inact %>%
  dplyr::filter(dose_frac == 2)



pool_log_ratio_inact_frac_0.5_day0 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 0)
pool_log_ratio_inact_frac_0.5_day14 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 14)
pool_log_ratio_inact_frac_0.5_day28 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 28)
pool_log_ratio_inact_frac_0.5_day32 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 32)
pool_log_ratio_inact_frac_0.5_day42 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 42)
pool_log_ratio_inact_frac_0.5_day56 = pool_log_ratio_frac_day(inact_frac_greater_0.5, 56)


pool_log_ratio_inact_full_day0 = pool_log_ratio_frac_day(inact_full, 0)
pool_log_ratio_inact_full_day14 = pool_log_ratio_frac_day(inact_full, 14)
pool_log_ratio_inact_full_day28 = pool_log_ratio_frac_day(inact_full, 28)
pool_log_ratio_inact_full_day32 = pool_log_ratio_frac_day(inact_full, 32)
pool_log_ratio_inact_full_day42 = pool_log_ratio_frac_day(inact_full, 42)
pool_log_ratio_inact_full_day56 = pool_log_ratio_frac_day(inact_full, 56)


pool_log_ratio_inact_frac_2_day0 = pool_log_ratio_frac_day(inact_frac_2, 0)
pool_log_ratio_inact_frac_2_day14 = pool_log_ratio_frac_day(inact_frac_2, 14)
pool_log_ratio_inact_frac_2_day28 = pool_log_ratio_frac_day(inact_frac_2, 28)
pool_log_ratio_inact_frac_2_day32 = pool_log_ratio_frac_day(inact_frac_2, 32)
pool_log_ratio_inact_frac_2_day42 = pool_log_ratio_frac_day(inact_frac_2, 42)
pool_log_ratio_inact_frac_2_day56 = pool_log_ratio_frac_day(inact_frac_2, 56)



## Protein subnit, fraction = 0.1, 0.3, 0.33, 0.5, 1, 1.5, 1.67, 2, 2.5, 5

data_protein = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="Protein subunit")

write.csv(data_protein, file = file.path(data_dir, "data_protein.csv"))

protein_frac_less_0.5 = data_protein %>%
  dplyr::filter(dose_frac < 0.5)

protein_frac_greater_0.5 = data_protein %>%
  dplyr::filter(dose_frac >= 0.5 & dose_frac < 1)

protein_full = data_protein %>%
  dplyr::filter(dose_frac == 1)

protein_frac_2 = data_protein %>%
  dplyr::filter(dose_frac == 2)


pool_log_ratio_protein_frac_less_0.5_day0 = pool_log_ratio_frac_day(protein_frac_less_0.5, 0)
pool_log_ratio_protein_frac_less_0.5_day14 = pool_log_ratio_frac_day(protein_frac_less_0.5, 14)
pool_log_ratio_protein_frac_less_0.5_day21 = pool_log_ratio_frac_day(protein_frac_less_0.5, 21)
pool_log_ratio_protein_frac_less_0.5_day28 = pool_log_ratio_frac_day(protein_frac_less_0.5, 28)
pool_log_ratio_protein_frac_less_0.5_day35 = pool_log_ratio_frac_day(protein_frac_less_0.5, 35)
pool_log_ratio_protein_frac_less_0.5_day42 = pool_log_ratio_frac_day(protein_frac_less_0.5, 42)
pool_log_ratio_protein_frac_less_0.5_day56 = pool_log_ratio_frac_day(protein_frac_less_0.5, 56)


pool_log_ratio_protein_frac_0.5_day21 = pool_log_ratio_frac_day(protein_frac_greater_0.5, 21)
pool_log_ratio_protein_frac_0.5_day28 = pool_log_ratio_frac_day(protein_frac_greater_0.5, 28)
pool_log_ratio_protein_frac_0.5_day35 = pool_log_ratio_frac_day(protein_frac_greater_0.5, 35)


pool_log_ratio_protein_full_day0 = pool_log_ratio_frac_day(protein_full, 0)
pool_log_ratio_protein_full_day14 = pool_log_ratio_frac_day(protein_full, 14)
pool_log_ratio_protein_full_day21 = pool_log_ratio_frac_day(protein_full, 21)
pool_log_ratio_protein_full_day28 = pool_log_ratio_frac_day(protein_full, 28)
pool_log_ratio_protein_full_day35 = pool_log_ratio_frac_day(protein_full, 35)
pool_log_ratio_protein_full_day42 = pool_log_ratio_frac_day(protein_full, 42)
pool_log_ratio_protein_full_day49 = pool_log_ratio_frac_day(protein_full, 49)
pool_log_ratio_protein_full_day56 = pool_log_ratio_frac_day(protein_full, 56)


pool_log_ratio_protein_frac_2_day0 = pool_log_ratio_frac_day(protein_frac_2, 0)
pool_log_ratio_protein_frac_2_day21 = pool_log_ratio_frac_day(protein_frac_2, 21)
pool_log_ratio_protein_frac_2_day35 = pool_log_ratio_frac_day(protein_frac_2, 35)
pool_log_ratio_protein_frac_2_day60 = pool_log_ratio_frac_day(protein_frac_2, 60)



## Non-replicating viral vector, fraction = 1, 2

data_viralvector = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="Non-replicating viral vector")

write.csv(data_viralvector, file = file.path(data_dir, "data_viralvector.csv"))

viralvector_full = data_viralvector %>%
  dplyr::filter(dose_frac == 1)

viralvector_frac_2 = data_viralvector %>%
  dplyr::filter(dose_frac == 2)



pool_log_ratio_viralvector_full_day0 = pool_log_ratio_frac_day(viralvector_full, 0)
pool_log_ratio_viralvector_full_day28 = pool_log_ratio_frac_day(viralvector_full, 28)
pool_log_ratio_viralvector_full_day56 = pool_log_ratio_frac_day(viralvector_full, 56)
pool_log_ratio_viralvector_full_day70 = pool_log_ratio_frac_day(viralvector_full, 70)


pool_log_ratio_viralvector_frac_2_day0 = pool_log_ratio_frac_day(viralvector_frac_2, 0)
pool_log_ratio_viralvector_frac_2_day28 = pool_log_ratio_frac_day(viralvector_frac_2, 28)
pool_log_ratio_viralvector_frac_2_day56 = pool_log_ratio_frac_day(viralvector_frac_2, 56)
pool_log_ratio_viralvector_frac_2_day70 = pool_log_ratio_frac_day(viralvector_frac_2, 70)



## Virus-like particle, fraction = 1, 2, 4

data_viralparticle = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="Virus-like particle")

write.csv(data_viralparticle, file = file.path(data_dir, "data_viralparticle.csv"))

viralparticle_full = data_viralparticle %>%
  dplyr::filter(dose_frac == 1)

viralparticle_frac_2 = data_viralparticle %>%
  dplyr::filter(dose_frac == 2)

viralparticle_frac_4 = data_viralparticle %>%
  dplyr::filter(dose_frac == 4)


pool_log_ratio_viralparticle_full_day0 = pool_log_ratio_frac_day(viralparticle_full, 0)
pool_log_ratio_viralparticle_full_day21 = pool_log_ratio_frac_day(viralparticle_full, 21)
pool_log_ratio_viralparticle_full_day42 = pool_log_ratio_frac_day(viralparticle_full, 42)


pool_log_ratio_viralparticle_frac_2_day0 = pool_log_ratio_frac_day(viralparticle_frac_2, 0)
pool_log_ratio_viralparticle_frac_2_day21 = pool_log_ratio_frac_day(viralparticle_frac_2, 21)
pool_log_ratio_viralparticle_frac_2_day42 = pool_log_ratio_frac_day(viralparticle_frac_2, 42)


pool_log_ratio_viralparticle_frac4_day0 = pool_log_ratio_frac_day(viralparticle_frac_4, 0)
pool_log_ratio_viralparticle_frac4_day21 = pool_log_ratio_frac_day(viralparticle_frac_4, 21)
pool_log_ratio_viralparticle_frac4_day42 = pool_log_ratio_frac_day(viralparticle_frac_4, 42)



## DNA based vaccine, fraction = 0.5, 1, 2 

data_DNA = data_trans %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(dose_frac) 
                & schedule_type == "prime" & vaccine_type =="DNA based vaccine")

write.csv(data_DNA, file = file.path(data_dir, "data_DNA.csv"))

DNA_frac_great_0.5 = data_DNA %>%
  dplyr::filter(dose_frac >= 0.5 & dose_frac < 1)

DNA_full = data_DNA %>%
  dplyr::filter(dose_frac == 1)

DNA_frac_2 = data_DNA %>%
  dplyr::filter(dose_frac == 2)


pool_log_ratio_DNA_frac_0.5_day0 = pool_log_ratio_frac_day(DNA_frac_great_0.5, 0)
pool_log_ratio_DNA_frac_0.5_day28 = pool_log_ratio_frac_day(DNA_frac_great_0.5, 28)
pool_log_ratio_DNA_frac_0.5_day56 = pool_log_ratio_frac_day(DNA_frac_great_0.5, 56)
pool_log_ratio_DNA_frac_0.5_day84 = pool_log_ratio_frac_day(DNA_frac_great_0.5, 84)


pool_log_ratio_DNA_full_day0 = pool_log_ratio_frac_day(DNA_full, 0)
pool_log_ratio_DNA_full_day28 = pool_log_ratio_frac_day(DNA_full, 28)
pool_log_ratio_DNA_full_day56 = pool_log_ratio_frac_day(DNA_full, 56)
pool_log_ratio_DNA_full_day84 = pool_log_ratio_frac_day(DNA_full, 84)


pool_log_ratio_DNA_frac_2_day0 = pool_log_ratio_frac_day(DNA_frac_2, 0)
pool_log_ratio_DNA_frac_2_day42 = pool_log_ratio_frac_day(DNA_frac_2, 42)
pool_log_ratio_DNA_frac_2_day56 = pool_log_ratio_frac_day(DNA_frac_2, 56)
pool_log_ratio_DNA_frac_2_day84 = pool_log_ratio_frac_day(DNA_frac_2, 84)
