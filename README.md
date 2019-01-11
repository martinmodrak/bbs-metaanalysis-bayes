# Bayesian metaanalysis of phenotype-genotype relationship in Bardet-Biedl Syndrome

This repository contains code and data accompanying the upcoming publication "Exploring the function of the BBSome using clinical data: Meta-analysis of genotype-phenotype associations in Bardet-Biedl Syndrome"

Contents of the files:

- `data` folder, contains the main dataset and results of the frequentist analyses from the main manuscript
- `main_analysis.Rmd` source for Appendix 2, describing (both in accessible and a mathy way) the model reported in the main manuscript. Reproduces all of the Bayesian figures from the main manscript + some additional checks and insights
- `alternative_models.Rmd` source for Appendix 4, describes all alternative models we tried throughout the analysis and shows posterior predictive checks that guided our selection of the model for the main analysis
- `multiverse_analysis.Rmd` source for Appendix 3, compares how the main conclusions of the paper hold under all alternative models, both Bayesian and frequentist.
- `data_processing.R` code to load and preprocess the dataset
- `models.R` definitions of all models used in the analysis
- `models_funcs.R` helper functions to easily fit the dataset with all models
- `plots.R` code for all the fancy plots used in the analysis