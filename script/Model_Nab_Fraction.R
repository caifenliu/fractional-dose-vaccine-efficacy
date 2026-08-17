library(readxl)
library(estmeansd)
library(tidyverse)
library(nlme)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(knitr)
library(scales)
library(lemon)
library(ggsci)
library(meta)



# Run from the `script` folder.
dir.create("Output", showWarnings = FALSE, recursive = TRUE)
neutralisation_data_dir <- file.path("..", "processed data", "dose neutralisation")
symptomatic_data_dir <- file.path("..", "processed data", "symptomatic efficacy")
asymptomatic_data_dir <- file.path("..", "processed data", "asymptomatic efficacy")

neutralisation_se_lookup <- read_csv(file.path(neutralisation_data_dir, "data_trans.csv"), show_col_types = FALSE)

# Gauss-Hermite quadrature for averaging protection over individual neutralisation levels.
n_quadrature_points <- 61
quadrature_index <- seq_len(n_quadrature_points - 1)
quadrature_matrix <- matrix(0, n_quadrature_points, n_quadrature_points)
quadrature_off_diagonal <- sqrt(quadrature_index / 2)
quadrature_matrix[cbind(quadrature_index, quadrature_index + 1)] <- quadrature_off_diagonal
quadrature_matrix[cbind(quadrature_index + 1, quadrature_index)] <- quadrature_off_diagonal
quadrature_eigen <- eigen(quadrature_matrix, symmetric = TRUE)
quadrature_order <- order(quadrature_eigen$values)
normal_quadrature_nodes <- sqrt(2) * quadrature_eigen$values[quadrature_order]
normal_quadrature_weights <- quadrature_eigen$vectors[1, quadrature_order]^2

logistic_protection_integrated <- function(mu, response_sd, k, c50) {
  out <- rep(NA, length(mu))
  ok <- !is.na(mu) & !is.na(response_sd) & response_sd >= 0 & !is.na(k) & !is.na(c50)
  if (!any(ok)) {
    return(out)
  }

  eta <- outer(mu[ok] - c50, response_sd * normal_quadrature_nodes, "+")
  out[ok] <- as.vector(plogis(k * eta) %*% normal_quadrature_weights)

  out
}

sym_efficacy_input <- read_csv(file.path(symptomatic_data_dir, "SummaryTable_Efficacy_NeutRatio_SD_SEM_Sym.csv"), show_col_types = FALSE)
asym_efficacy_input <- read_csv(file.path(asymptomatic_data_dir, "Summary_Efficacy_NeutRatio_SD_SEM_Asym.csv"), show_col_types = FALSE)
asym_efficacy_input$SeSD[3] <- mean((asym_efficacy_input %>% dplyr::filter(!is.na(SeSD)))$SeSD)

sym_response_sd <- sym_efficacy_input$PooledSD[!is.na(sym_efficacy_input$PooledSD)][1]
asym_response_sd <- metagen(
  TE = SD,
  seTE = SeSD,
  studlab = Study,
  data = asym_efficacy_input,
  n.e = NumberIndividuals_Vaccine,
  level.ci = 0.95,
  hakn = TRUE
)$TE.random

sym_k <- 3.09770274348263
sym_c50 <- log10(0.201088532847589)
asym_k <- 1.83949974964737
asym_c50 <- log10(2.6223245183247)

add_integrated_efficacy <- function(data) {
  data$Efficacy_asym <- logistic_protection_integrated(
    data$Log_ratio_to_conv,
    asym_response_sd,
    asym_k,
    asym_c50
  )
  data$Efficacy_asym_lower <- logistic_protection_integrated(
    data$Log_ratio_to_conv_lower,
    asym_response_sd,
    asym_k,
    asym_c50
  )
  data$Efficacy_asym_upper <- logistic_protection_integrated(
    data$Log_ratio_to_conv_upper,
    asym_response_sd,
    asym_k,
    asym_c50
  )

  data$Efficacy_sym <- logistic_protection_integrated(
    data$Log_ratio_to_conv,
    sym_response_sd,
    sym_k,
    sym_c50
  )
  data$Efficacy_sym_lower <- logistic_protection_integrated(
    data$Log_ratio_to_conv_lower,
    sym_response_sd,
    sym_k,
    sym_c50
  )
  data$Efficacy_sym_upper <- logistic_protection_integrated(
    data$Log_ratio_to_conv_upper,
    sym_response_sd,
    sym_k,
    sym_c50
  )

  data
}

attach_logratio_se <- function(data) {
  if ("se_log_ratio_to_conv" %in% names(data) && all(!is.na(data$se_log_ratio_to_conv))) {
    return(data)
  }

  join_keys <- intersect(
    c("study_id", "study", "vaccine", "vaccine_type", "dose_frac", "measured_day", "N", "nab", "conv"),
    intersect(names(data), names(neutralisation_se_lookup))
  )

  data$.join_key <- do.call(paste, c(lapply(data[join_keys], as.character), sep = "\r"))
  lookup_data <- neutralisation_se_lookup
  lookup_data$.join_key <- do.call(paste, c(lapply(lookup_data[join_keys], as.character), sep = "\r"))

  lookup <- lookup_data %>%
    dplyr::select(
      .join_key,
      se_log_ratio_to_conv,
      se_conv,
      sd_conv
    ) %>%
    dplyr::group_by(.join_key) %>%
    dplyr::summarise(
      se_log_ratio_to_conv = se_log_ratio_to_conv[!is.na(se_log_ratio_to_conv)][1],
      se_conv = se_conv[!is.na(se_conv)][1],
      sd_conv = sd_conv[!is.na(sd_conv)][1],
      .groups = "drop"
    )

  data %>%
    dplyr::select(-dplyr::any_of(c("se_log_ratio_to_conv", "se_conv", "sd_conv"))) %>%
    dplyr::left_join(lookup, by = ".join_key") %>%
    dplyr::select(-.join_key)
}

### Load peak neutralisation data for dose-response models
mRNA_peak = read_csv(file.path(neutralisation_data_dir, "data_peak_mRNA.csv")) %>% attach_logratio_se()
mRNA_peak = mRNA_peak %>% dplyr::filter(dose_frac <= 1 & dose_frac != 0.5)

inact_peak = read_csv(file.path(neutralisation_data_dir, "data_peak_inact.csv")) %>% attach_logratio_se()

protein_peak = read_csv(file.path(neutralisation_data_dir, "data_peak_protein.csv")) %>% attach_logratio_se()
protein_peak = protein_peak %>% dplyr::filter(dose_frac < 2) 

viralvector_peak = read_csv(file.path(neutralisation_data_dir, "data_peak_viralvector.csv")) %>% attach_logratio_se()
viralvector_peak <- viralvector_peak %>%
  dplyr::mutate(
    log_ratio_se = se_log_ratio_to_conv,
    log_ratio_lower = log_ratio_to_conv - 1.96 * log_ratio_se,
    log_ratio_upper = log_ratio_to_conv + 1.96 * log_ratio_se
  ) %>%
  dplyr::filter(!is.na(log_ratio_to_conv) & !is.na(log_ratio_se))

## Fit weighted dose-response models with shared convalescent covariance.
##
## The log-ratio SE is the marginal uncertainty of each study arm:
## sqrt(SE_vaccine^2 + SE_convalescent^2). For estimating the dose-response
## slope, vaccine arms that share the same convalescent comparator should also
## share the convalescent error term. The GLS covariance matrix therefore has
## independent vaccine-arm uncertainty on the diagonal and a shared convalescent
## covariance block for rows with the same convalescent comparator.

frac_pred <- data.frame("dose_frac" = seq(0.03, 5, by = 0.001))

make_convalescent_group <- function(data) {
  group_keys <- intersect(
    c("study_id", "study", "ClinicalTrials", "against_strain", "assay", "conv", "conv_lower", "conv_upper", "conv_N"),
    names(data)
  )
  do.call(paste, c(lapply(data[group_keys], as.character), sep = "\r"))
}

prepare_gls_data <- function(data) {
  data %>%
    dplyr::mutate(
      log_ratio_se = dplyr::coalesce(
        dplyr::if_else(!is.na(se_log_ratio_to_conv), se_log_ratio_to_conv, NA),
        sqrt(se_vaccine^2 + se_conv^2)
      ),
      convalescent_group = make_convalescent_group(.)
    ) %>%
    dplyr::filter(
      !is.na(log_ratio_to_conv),
      !is.na(dose_frac),
      dose_frac > 0,
      !is.na(se_vaccine),
      !is.na(se_conv),
      se_vaccine > 0,
      se_conv >= 0
    )
}

build_shared_convalescent_covariance <- function(data) {
  n <- nrow(data)
  sigma <- diag(data$se_vaccine^2, nrow = n, ncol = n)

  for (group in unique(data$convalescent_group)) {
    idx <- which(data$convalescent_group == group)
    conv_var <- data$se_conv[idx][!is.na(data$se_conv[idx])][1]^2
    if (!is.na(conv_var) && conv_var > 0) {
      sigma[idx, idx] <- sigma[idx, idx] + conv_var
    }
  }

  sigma + diag(1e-10, nrow = n, ncol = n)
}

fit_weighted_dose_response <- function(data, platform_label) {
  model_data <- prepare_gls_data(data)
  if (nrow(model_data) < 3) {
    stop(platform_label, " has fewer than 3 observations available for GLS fitting.")
  }

  x <- log10(model_data$dose_frac)
  y <- model_data$log_ratio_to_conv
  design <- cbind("(Intercept)" = 1, "log10(dose_frac)" = x)
  covariance <- build_shared_convalescent_covariance(model_data)

  sigma_inv_design <- solve(covariance, design)
  sigma_inv_y <- solve(covariance, y)
  information <- t(design) %*% sigma_inv_design
  beta_covariance <- solve(information)
  beta <- as.vector(beta_covariance %*% t(design) %*% sigma_inv_y)
  names(beta) <- colnames(design)

  fitted <- as.vector(design %*% beta)
  residual <- y - fitted
  weighted_sse <- as.numeric(t(residual) %*% solve(covariance, residual))
  df <- nrow(model_data) - ncol(design)
  dispersion <- ifelse(df > 0, max(1, weighted_sse / df), 1)
  beta_covariance <- beta_covariance * dispersion

  weighted_mean <- sum(y / diag(covariance)) / sum(1 / diag(covariance))
  weighted_sst <- sum((y - weighted_mean)^2 / diag(covariance))
  weighted_r2 <- ifelse(weighted_sst > 0, 1 - sum((y - fitted)^2 / diag(covariance)) / weighted_sst, NA)
  covariance_inv <- solve(covariance)
  intercept_design <- matrix(1, nrow = length(y), ncol = 1)
  intercept_only <- as.numeric(
    solve(t(intercept_design) %*% covariance_inv %*% intercept_design) %*%
      t(intercept_design) %*% covariance_inv %*% y
  )
  intercept_residual <- y - intercept_only
  gls_sst <- as.numeric(t(intercept_residual) %*% covariance_inv %*% intercept_residual)
  gls_r2 <- ifelse(gls_sst > 0, 1 - weighted_sse / gls_sst, NA)
  slope_p_value <- ifelse(
    df > 0,
    2 * stats::pt(abs(beta[2] / sqrt(beta_covariance[2, 2])), df = df, lower.tail = FALSE),
    NA
  )

  list(
    platform = platform_label,
    data = model_data,
    beta = beta,
    beta_covariance = beta_covariance,
    dispersion = dispersion,
    r_squared = weighted_r2,
    gls_r_squared = gls_r2,
    slope_p_value = slope_p_value,
    n = nrow(model_data),
    n_convalescent_groups = dplyr::n_distinct(model_data$convalescent_group)
  )
}

format_gls_label <- function(fit) {
  p_text <- dplyr::case_when(
    is.na(fit$slope_p_value) ~ "P==NA",
    fit$slope_p_value < 0.001 ~ "P<0.001",
    TRUE ~ paste0("P==", sprintf("%.3f", fit$slope_p_value))
  )
  paste0("atop(R^2==", sprintf("%.3f", fit$gls_r_squared), ",", p_text, ")")
}

predict_weighted_dose_response <- function(fit, new_data = frac_pred) {
  design <- cbind("(Intercept)" = 1, "log10(dose_frac)" = log10(new_data$dose_frac))
  estimate <- as.vector(design %*% fit$beta)
  se <- sqrt(rowSums((design %*% fit$beta_covariance) * design))
  data.frame(
    Dose_fraction = new_data$dose_frac,
    Log_ratio_to_conv = estimate,
    Log_ratio_to_conv_lower = estimate - 1.96 * se,
    Log_ratio_to_conv_upper = estimate + 1.96 * se
  )
}

fit_weighted_mRNA <- fit_weighted_dose_response(mRNA_peak, "mRNA")
fit_weighted_inact <- fit_weighted_dose_response(inact_peak, "Inactivated vaccine")
fit_weighted_protein <- fit_weighted_dose_response(protein_peak, "Protein subunit")
fit_weighted_viralvector <- fit_weighted_dose_response(viralvector_peak, "Viral vector")

make_raw_logratio_plot_data <- function(fit) {
  fit$data %>%
    dplyr::mutate(
      log_ratio_lower = log_ratio_to_conv - 1.96 * log_ratio_se,
      log_ratio_upper = log_ratio_to_conv + 1.96 * log_ratio_se
    )
}

mRNA_logratio_points <- make_raw_logratio_plot_data(fit_weighted_mRNA)
inact_logratio_points <- make_raw_logratio_plot_data(fit_weighted_inact)
protein_logratio_points <- make_raw_logratio_plot_data(fit_weighted_protein)
viralvector_logratio_points <- make_raw_logratio_plot_data(fit_weighted_viralvector)

fraction_nab_fit_summary <- dplyr::bind_rows(
  data.frame(
    vaccine_type = "mRNA",
    method = "weighted_gls_shared_convalescent_covariance",
    n = fit_weighted_mRNA$n,
    n_convalescent_groups = fit_weighted_mRNA$n_convalescent_groups,
    intercept = fit_weighted_mRNA$beta[1],
    intercept_se = sqrt(fit_weighted_mRNA$beta_covariance[1, 1]),
    slope = fit_weighted_mRNA$beta[2],
    slope_se = sqrt(fit_weighted_mRNA$beta_covariance[2, 2]),
    weighted_r_squared = fit_weighted_mRNA$r_squared,
    gls_r_squared = fit_weighted_mRNA$gls_r_squared,
    slope_p_value = fit_weighted_mRNA$slope_p_value,
    dispersion = fit_weighted_mRNA$dispersion
  ),
  data.frame(
    vaccine_type = "Inactivated vaccine",
    method = "weighted_gls_shared_convalescent_covariance",
    n = fit_weighted_inact$n,
    n_convalescent_groups = fit_weighted_inact$n_convalescent_groups,
    intercept = fit_weighted_inact$beta[1],
    intercept_se = sqrt(fit_weighted_inact$beta_covariance[1, 1]),
    slope = fit_weighted_inact$beta[2],
    slope_se = sqrt(fit_weighted_inact$beta_covariance[2, 2]),
    weighted_r_squared = fit_weighted_inact$r_squared,
    gls_r_squared = fit_weighted_inact$gls_r_squared,
    slope_p_value = fit_weighted_inact$slope_p_value,
    dispersion = fit_weighted_inact$dispersion
  ),
  data.frame(
    vaccine_type = "Protein subunit",
    method = "weighted_gls_shared_convalescent_covariance",
    n = fit_weighted_protein$n,
    n_convalescent_groups = fit_weighted_protein$n_convalescent_groups,
    intercept = fit_weighted_protein$beta[1],
    intercept_se = sqrt(fit_weighted_protein$beta_covariance[1, 1]),
    slope = fit_weighted_protein$beta[2],
    slope_se = sqrt(fit_weighted_protein$beta_covariance[2, 2]),
    weighted_r_squared = fit_weighted_protein$r_squared,
    gls_r_squared = fit_weighted_protein$gls_r_squared,
    slope_p_value = fit_weighted_protein$slope_p_value,
    dispersion = fit_weighted_protein$dispersion
  ),
  data.frame(
    vaccine_type = "Viral vector",
    method = "weighted_gls_shared_convalescent_covariance",
    n = fit_weighted_viralvector$n,
    n_convalescent_groups = fit_weighted_viralvector$n_convalescent_groups,
    intercept = fit_weighted_viralvector$beta[1],
    intercept_se = sqrt(fit_weighted_viralvector$beta_covariance[1, 1]),
    slope = fit_weighted_viralvector$beta[2],
    slope_se = sqrt(fit_weighted_viralvector$beta_covariance[2, 2]),
    weighted_r_squared = fit_weighted_viralvector$r_squared,
    gls_r_squared = fit_weighted_viralvector$gls_r_squared,
    slope_p_value = fit_weighted_viralvector$slope_p_value,
    dispersion = fit_weighted_viralvector$dispersion
  )
)

write.csv(fraction_nab_fit_summary, "Output/fraction_nab_fit_summary.csv", row.names = FALSE)


### Get weighted GLS predictions for inactivated vaccine diagnostic plot

Model_frac_nab_inactivated <- predict_weighted_dose_response(fit_weighted_inact)
Model_frac_nab_inactivated <- add_integrated_efficacy(Model_frac_nab_inactivated)
write.csv(Model_frac_nab_inactivated, "Output/model_frac_nab_inactivated.csv", row.names = FALSE)


## Plot fraction-logratio for inactivated vaccine
Fit_frac_logratio_inactivated <- ggplot()+
  geom_ribbon(data = Model_frac_nab_inactivated, aes(ymin = Log_ratio_to_conv_lower, ymax = Log_ratio_to_conv_upper, x = Dose_fraction), fill = "#F4A582") +
  geom_line(data = Model_frac_nab_inactivated, aes(x = Dose_fraction, y = Log_ratio_to_conv), colour = "#D73027") +
  geom_errorbar(data = inact_logratio_points, aes(x = dose_frac, ymin = log_ratio_lower, ymax = log_ratio_upper),
                width = 0.025, linewidth = 0.25, colour = "grey35", alpha = 0.35) +
  geom_point(data = inact_logratio_points, aes(x = dose_frac, y = log_ratio_to_conv),
             colour = "#D73027", alpha = 0.55, size = 1.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(limits = c(0,2)) +
  scale_y_continuous(limits = c(-1.75, 1.5)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  annotate("text",x=1.8,y=1.3,label=format_gls_label(fit_weighted_inact),parse=T, size = 4.5) + labs(title = "Inactivated vaccine", x = "Dose fraction", y = "Log normalised neutralisation level") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Log normalised neutralisation level") + theme(axis.title.y = element_text(size = 15))

## Plot fraction-efficacy for inactivated vaccine
Fit_frac_efficacy_inactivated_sym<-ggplot(data=Model_frac_nab_inactivated, aes(y=100*Efficacy_sym,x=Dose_fraction)) +
  geom_ribbon(data=Model_frac_nab_inactivated,aes(ymin=100*Efficacy_sym_lower, ymax=100*Efficacy_sym_upper), fill = "#F4A582")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_inactivated,aes(y=100*Efficacy_sym,x=Dose_fraction), colour = "#D73027") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against symptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))

Fit_frac_efficacy_inactivated_asym<-ggplot(data=Model_frac_nab_inactivated, aes(y=100*Efficacy_asym,x=Dose_fraction)) +
  geom_ribbon(data=Model_frac_nab_inactivated,aes(ymin=100*Efficacy_asym_lower, ymax=100*Efficacy_asym_upper), fill = "#F4A582")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), labels=c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_inactivated,aes(y=100*Efficacy_asym,x=Dose_fraction), colour = "#D73027") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against asymptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))

### Get weighted GLS predictions for mRNA

Model_frac_nab_mRNA <- predict_weighted_dose_response(fit_weighted_mRNA)
Model_frac_nab_mRNA <- add_integrated_efficacy(Model_frac_nab_mRNA)
write.csv(Model_frac_nab_mRNA, "Output/model_frac_nab_mRNA.csv", row.names = FALSE)


## Plot fraction-logratio for mRNA
Fit_frac_logratio_mRNA <- ggplot()+
  geom_ribbon(data = Model_frac_nab_mRNA, aes(ymin = Log_ratio_to_conv_lower, ymax = Log_ratio_to_conv_upper, x = Dose_fraction), fill = "#F2C14E") +
  geom_line(data = Model_frac_nab_mRNA, aes(x = Dose_fraction, y = Log_ratio_to_conv), colour = "#8C510A") +
  geom_errorbar(data = mRNA_logratio_points, aes(x = dose_frac, ymin = log_ratio_lower, ymax = log_ratio_upper),
                width = 0.025, linewidth = 0.25, colour = "grey35", alpha = 0.35) +
  geom_point(data = mRNA_logratio_points, aes(x = dose_frac, y = log_ratio_to_conv),
             colour = "#8C510A", alpha = 0.55, size = 1.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(limits = c(0,2)) +
  scale_y_continuous(limits = c(-1.75, 1.5)) + 
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  annotate("text",x=1.8,y=1.3,label=format_gls_label(fit_weighted_mRNA),parse=T, size = 4.5) + labs(title = "mRNA", x = "Dose fraction", y = "Log normalised neutralisation level") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Log normalised neutralisation level") + theme(axis.title.y = element_text(size = 15))  

## Plot fraction-efficacy (sym, asym, hospital, severe) for mRNA
Fit_frac_efficacy_mRNA_sym<-ggplot(data=Model_frac_nab_mRNA, aes(y=100*Efficacy_sym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_mRNA,aes(ymin=100*Efficacy_sym_lower, ymax=100*Efficacy_sym_upper), fill = "#F2C14E")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_mRNA,aes(y=100*Efficacy_sym,x=Dose_fraction), colour = "#8C510A") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against symptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))

Fit_frac_efficacy_mRNA_asym<-ggplot(data=Model_frac_nab_mRNA, aes(y=100*Efficacy_asym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_mRNA,aes(ymin=100*Efficacy_asym_lower, ymax=100*Efficacy_asym_upper), fill = "#F2C14E")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_mRNA,aes(y=100*Efficacy_asym,x=Dose_fraction), colour = "#8C510A") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against asymptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))
  

### Get weighted GLS predictions for protein subunit

Model_frac_nab_protein <- predict_weighted_dose_response(fit_weighted_protein)
Model_frac_nab_protein <- add_integrated_efficacy(Model_frac_nab_protein)
write.csv(Model_frac_nab_protein, "Output/model_frac_nab_protein.csv", row.names = FALSE)


## Plot fraction-logratio for protein subunit

Fit_frac_logratio_protein <- ggplot()+
  geom_ribbon(data = Model_frac_nab_protein, aes(ymin = Log_ratio_to_conv_lower, ymax = Log_ratio_to_conv_upper, x = Dose_fraction), fill = "#C7EAE5") +
  geom_line(data = Model_frac_nab_protein, aes(x = Dose_fraction, y = Log_ratio_to_conv), colour = "#01665E") +
  geom_errorbar(data = protein_logratio_points, aes(x = dose_frac, ymin = log_ratio_lower, ymax = log_ratio_upper),
                width = 0.025, linewidth = 0.25, colour = "grey35", alpha = 0.35) +
  geom_point(data = protein_logratio_points, aes(x = dose_frac, y = log_ratio_to_conv),
             colour = "#01665E", alpha = 0.55, size = 1.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(-1.75, 1.5)) + 
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  annotate("text",x=1.8,y=1.3,label=format_gls_label(fit_weighted_protein),parse=T, size = 4.5) + labs(title = "Protein subunit") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Log normalised neutralisation level") + theme(axis.title.y = element_text(size = 15))  

## Plot fraction-efficacy (sym, asym, hospital, severe) for protein subunit

Fit_frac_efficacy_protein_sym<-ggplot(data=Model_frac_nab_protein, aes(y=100*Efficacy_sym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_protein,aes(ymin=100*Efficacy_sym_lower, ymax=100*Efficacy_sym_upper), fill = "#C7EAE5")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_protein,aes(y=100*Efficacy_sym,x=Dose_fraction), colour = "#01665E") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against symptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))


Fit_frac_efficacy_protein_asym<-ggplot(data=Model_frac_nab_protein, aes(y=100*Efficacy_asym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_protein,aes(ymin=100*Efficacy_asym_lower, ymax=100*Efficacy_asym_upper), fill = "#C7EAE5")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_protein,aes(y=100*Efficacy_asym,x=Dose_fraction), colour = "#01665E") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against asymptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))
  

### Get weighted GLS predictions for viral vector

Model_frac_nab_viralvector <- predict_weighted_dose_response(fit_weighted_viralvector)
Model_frac_nab_viralvector <- add_integrated_efficacy(Model_frac_nab_viralvector)
write.csv(Model_frac_nab_viralvector, "Output/model_frac_nab_viralvector.csv", row.names = FALSE)


## Plot fraction-logratio for viral vector

Fit_frac_logratio_viralvector <- ggplot()+
  geom_ribbon(data = Model_frac_nab_viralvector, aes(ymin = Log_ratio_to_conv_lower, ymax = Log_ratio_to_conv_upper, x = Dose_fraction), fill = "#D8BFD8") +
  geom_line(data = Model_frac_nab_viralvector, aes(x = Dose_fraction, y = Log_ratio_to_conv), colour = "#5E3C99") +
  geom_errorbar(data = viralvector_logratio_points, aes(x = dose_frac, ymin = log_ratio_lower, ymax = log_ratio_upper),
                width = 0.025, linewidth = 0.25, colour = "grey35", alpha = 0.35) +
  geom_point(data = viralvector_logratio_points, aes(x = dose_frac, y = log_ratio_to_conv),
             colour = "#5E3C99", alpha = 0.55, size = 1.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(limits = c(0,2)) +
  scale_y_continuous(limits = c(-1.75, 1.5)) + 
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  annotate("text",x=1.8,y=1.3,label=format_gls_label(fit_weighted_viralvector),parse=T, size = 4.5) + labs(title = "NR viral vector") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Log normalised neutralisation level") + theme(axis.title.y = element_text(size = 15))  

## Plot fraction-efficacy (sym, asym, hospital, severe) for viral vector

Fit_frac_efficacy_viralvector_sym<-ggplot(data=Model_frac_nab_viralvector, aes(y=100*Efficacy_sym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_viralvector,aes(ymin=100*Efficacy_sym_lower, ymax=100*Efficacy_sym_upper), fill = "#D8BFD8")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), label = c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_viralvector,aes(y=100*Efficacy_sym,x=Dose_fraction), colour = "#5E3C99") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against symptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))
  

Fit_frac_efficacy_viralvector_asym<-ggplot(data=Model_frac_nab_viralvector, aes(y=100*Efficacy_asym,x=Dose_fraction)) +
  #adding the bands
  geom_ribbon(data=Model_frac_nab_viralvector,aes(ymin=100*Efficacy_asym_lower, ymax=100*Efficacy_asym_upper), fill = "#D8BFD8")+
  scale_x_log10(lim=c(0.03,5),breaks=c(0.03, 0.1,0.5,1,5),labels=c(0.03, 0.1, 0.5, 1, 5)) +
  scale_y_continuous(lim=c(0,100),breaks = seq(0, 100, by = 20), labels=c(0, 20, 40, 60, 80, 100)) +
  theme_linedraw() +
  theme(panel.grid.major=element_line(colour= "grey 80"),
        panel.background = element_rect(fill = "white",colour = NA),
        plot.background = element_rect(fill = "white",colour = NA),
        panel.grid.minor = element_blank(), legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.text = element_text(size = 15, colour = "black")) +
  geom_line(data=Model_frac_nab_viralvector,aes(y=100*Efficacy_asym,x=Dose_fraction), colour = "#5E3C99") +
  xlab("Dose fraction") + theme(axis.title.x = element_text(size = 15)) +
  ylab("Efficacy against asymptomatic infection (%)") + theme(axis.title.y = element_text(size = 15))


### Layout plots

layout_frac_logratio_2x2 <- (
  (Fit_frac_logratio_mRNA + labs(title = "mRNA") |
     Fit_frac_logratio_protein + labs(title = "Protein subunit")) /
    (Fit_frac_logratio_viralvector + labs(title = "NR viral vector") |
       Fit_frac_logratio_inactivated + labs(title = "Inactivated vaccine"))
) &
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 13, colour = "black"),
    plot.margin = margin(t = 16, r = 5, b = 5, l = 18)
)

pdf("Output/Layout_frac_logratio_2x2.pdf", height = 9, width = 9.4)
print(layout_frac_logratio_2x2)
dev.off()
png("Output/Layout_frac_logratio_2x2.png", height = 9, width = 9.4, units = "in", res = 300)
print(layout_frac_logratio_2x2)
dev.off()


layout_frac_efficacy <- (
  (Fit_frac_efficacy_mRNA_sym + labs(title = "mRNA", x = NULL, tag = "A.") |
     Fit_frac_efficacy_protein_sym + labs(title = "Protein subunit", x = NULL) |
     Fit_frac_efficacy_viralvector_sym + labs(title = "NR viral vector", x = NULL) |
     Fit_frac_efficacy_inactivated_sym + labs(title = "Inactivated vaccine", x = NULL)) /
    (Fit_frac_efficacy_mRNA_asym + labs(title = NULL, tag = "B.") |
       Fit_frac_efficacy_protein_asym + labs(title = NULL) |
       Fit_frac_efficacy_viralvector_asym + labs(title = NULL) |
       Fit_frac_efficacy_inactivated_asym + labs(title = NULL))
) &
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12.5, colour = "black"),
    axis.title.y = element_text(size = 12.5, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.tag = element_text(size = 16, face = "bold", family = "Times", colour = "black"),
    plot.tag.position = c(-0.035, 1.06),
    plot.margin = margin(t = 16, r = 5, b = 5, l = 18)
  )

pdf("Output/Layout_frac_efficacy.pdf", height = 8, width = 16)
print(layout_frac_efficacy)
dev.off()
