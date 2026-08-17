# Projection of Fractional-Dose COVID-19 Vaccine Efficacy Using Neutralising Antibody Levels

This repository contains the processed analysis data and R code used for Chapter 8 of Caifen Liu's PhD thesis.

The analyses in this chapter build on the methods and findings reported in the following published work:

> Du Z, Liu C, Bai Y, Wang L, Lim WW, Lau EHY, Cowling BJ.  
> **Predicting Efficacies of Fractional Doses of Vaccines by Using Neutralizing Antibody Levels: Systematic Review and Meta-Analysis.**  
> *JMIR Public Health and Surveillance.* 2024;10:e49812.  
> https://doi.org/10.2196/49812

## Repository structure

### `processed data/`

- `dose neutralisation/`: Transformed, vaccine-platform-specific, and peak neutralisation datasets used in the dose–neutralisation models.
- `symptomatic efficacy/`: Analysis-ready data used in the symptomatic neutralisation–efficacy model and selected fractional-dose projections.
- `asymptomatic efficacy/`: Analysis-ready data used in the asymptomatic neutralisation–efficacy model and selected fractional-dose projections.

### `script/`

- `Frac_DataTrans.R`: Preprocesses the dose–neutralisation data and generates transformed datasets.
- `Model_Nab_Fraction.R`: Runs the dose–neutralisation models, projects fractional-dose vaccine efficacy, and produces the main plots.
- `Combine_full_fit_plots.R`: Combines the neutralisation–efficacy plots.
- `Symptomatic/Model_Efficacy_Nab_Symptomatic.R`: Fits the neutralisation–efficacy model for symptomatic infection.
- `Asymptomatic/Model_Efficacy_Nab_Asymptomatic.R`: Fits the neutralisation–efficacy model for asymptomatic infection.

## Running the analysis

Run the modelling scripts from the `script/` directory. The scripts read their inputs from `../processed data/`.

The recommended run order is:

1. `Symptomatic/Model_Efficacy_Nab_Symptomatic.R`
2. `Asymptomatic/Model_Efficacy_Nab_Asymptomatic.R`
3. `Model_Nab_Fraction.R`
4. `Combine_full_fit_plots.R`

## Citation

If you use this code or the associated processed data, please cite:

Du Z, Liu C, Bai Y, Wang L, Lim WW, Lau EHY, Cowling BJ. Predicting Efficacies of Fractional Doses of Vaccines by Using Neutralizing Antibody Levels: Systematic Review and Meta-Analysis. *JMIR Public Health and Surveillance.* 2024;10:e49812. https://doi.org/10.2196/49812
