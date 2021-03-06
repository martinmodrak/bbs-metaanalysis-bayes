# Bayesian metaanalysis of phenotype-genotype relationship in Bardet-Biedl Syndrome

This repository contains code and data accompanying the upcoming publication "Exploring the function of the BBSome using clinical data: Meta-analysis of genotype-phenotype associations in Bardet-Biedl Syndrome"

Contents of the files:

- `data` folder, contains the main dataset and results of the frequentist analyses from the main manuscript
- `main_analysis.Rmd` describes (both in accessible and a mathy way) the model reported in the main manuscript. Reproduces all of the Bayesian figures from the main manscript + some additional checks and insights
- `alternative_models.Rmd` describes all alternative models we tried throughout the analysis and shows posterior predictive checks that guided our selection of the model for the main analysis
- `multiverse_analysis.Rmd` compares how the main conclusions of the paper hold under all alternative models, both Bayesian and frequentist.
- `validation_dataset.Rmd` code made to validate the findings on new data. We were later denied access to new data and the code here is discontinued.
- `master_document.Rmd` includes the 3 Rmd files above to form a single document for journal submissions
- `master_document_with_validation.Rmd` an alternative to `mater_document.Rmd` also including validation. Discontinued after we were denied access to the new dataset.
- `data_processing.R` code to load and preprocess the dataset
- `models.R` definitions of all models used in the analysis
- `models_funcs.R` helper functions to easily fit the dataset with all models
- `plots.R` code for all the fancy plots used in the analysis

## Rerunning the analysis


To rerun, you need to install `brms` whic requires `rstan`. Installing `rstan` directly with `install.package` may fail on some systems - see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started for more details.

The model fitting for the analyses may be time consuming. To avoid recomputing, fits are stored in the `stored_fits` directory, which is created on the fly. To save you time refitting, you can download the fits we used at: https://zenodo.org/record/3243270 (DOI: 10.5281/zenodo.3243270). Downloading the fits should also let you rerun the code without completely configuring `rstan`.


You can rerun the complete analysis be knitting `master_document.Rmd`, but the individual parts can also be executed/knitted separately.

`main_analysis.Rmd` is self-sufficient and can be run directly. It fits the main model, which should be relatively quick  (20min - 1 hour).

`alternative_models.Rmd` is self-sufficient and can be run directly. This will fit a large number of relatively large models and may take upwards of a day even on a powerful machine.

To run `multiverse_analysis.Rmd` you need to run `alternative_models.Rmd` first (or download fitted models), but after that it is pretty quick.  
