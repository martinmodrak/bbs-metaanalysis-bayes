fit_base_model <- function(def, data_long) {
  if(is.null(def)) {
    stop("Invalid def")
  }
  
  if(def$imputed) {
    stop("Cannot be used on models with imputed data")
  }
  
  stored_fits_dir <- "stored_fits"
  if(!dir.exists(here(stored_fits_dir))) {
    dir.create(here(stored_fits_dir))
  }

  data <- switch (def$filter,
    "none" = data_long,
    "lof" = data_long %>% filter(loss_of_function == "certain"),
    "age" = data_long %>% filter(!is.na(age_std_for_model)),
    "sex" = data_long %>% filter(!is.na(Sex)),
    stop("Unrecognized filter")
  )
  
  brm(formula = def$formula, data = data, prior = def$priors, file = here(stored_fits_dir,def$name),
      control = list(adapt_delta = 0.95))   
}

data_for_prediction_base_model <- function(def, genes_to_show, phenotypes_to_show, age = 30, 
                                           sex = "F") {
  result <-  tibble(gene = genes_to_show) %>%
    crossing(tibble(phenotype = phenotypes_to_show)) %>% 
    mutate(functional_group = functional_group_for_gene(gene))
  
  if("source" %in% def$additional_components) {
    result$source = "new_source"
  }
  
  if("lof" %in% def$additional_components) {
    result$loss_of_function_certain = 1
  }
  
  if("age" %in% def$additional_components) {
    result <- result %>% crossing(tibble(age_std_for_model = age_transform(age)))
  }
  
  if("sex" %in% def$additional_components) {
    result <- result %>% crossing(tibble(Sex = sex))
  }
  
  result
}