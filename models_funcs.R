filter_data_by_model_def <- function(def, data_long) {
  switch (def$filter,
          "none" = data_long,
          "lof" = data_long %>% filter(loss_of_function == "certain"),
          "age" = data_long %>% filter(!is.na(age_std_for_model)),
          "sex" = data_long %>% filter(!is.na(Sex)),
          "age_sex" = data_long %>% filter(!is.na(age_std_for_model), !is.na(Sex)),
          stop("Unrecognized filter")
  )  
}

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

  data <- filter_data_by_model_def(def, data_long)
  
  brm(formula = def$formula, data = data, prior = def$priors, file = here(stored_fits_dir,def$name),
      control = list(adapt_delta = 0.95))   
}

fit_imputed_model <- function(def, data_imputed) {
  if(is.null(def)) {
    stop("Invalid def")
  }
  
  if(!def$imputed) {
    stop("Can only be used on models with imputed data")
  }
  
  stored_fits_dir <- "stored_fits"
  if(!dir.exists(here(stored_fits_dir))) {
    dir.create(here(stored_fits_dir))
  }
  
  brm_multiple(formula = def$formula, data = data_imputed, prior = def$priors, file = here(stored_fits_dir,def$name),
      control = list(adapt_delta = 0.95))   
}

data_for_prediction_base_model <- function(def, genes_to_show, phenotypes_to_show, ages = c(10,30,50), 
                                           sexes = c("F","M")) {
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
    result <- result %>% crossing(tibble(age_std_for_model = age_transform(ages)))
  }
  
  if("sex" %in% def$additional_components) {
    result <- result %>% crossing(tibble(Sex = sexes))
  }
  
  result
}

filter_for_BBSome <- function(data_for_prediction) {
  filter(data_for_prediction, functional_group == "BBSome", gene != "BBS18")
}

get_tidy_samples <- function(fit, data_for_prediction) {
  fitted_detailed <- fitted(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as.tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) %>%
    mutate(odds = value/(1-value))
  
}

get_samples_pairwise_diff <- function(model_def, samples_tidy) {
  component_to_join_map <- function(component) {
    switch(component,
           "lof" = c("loss_of_function_certain" = "loss_of_function_certain"),
           "age" = c("age_std_for_model" = "age_std_for_model"),
           "sex" = c("Sex" = "Sex"),
           "source" = c()
    )
  }
  
  join_additional <- model_def$additional_components %>% map(component_to_join_map) %>% do.call(c, .)
  join_by <- c("sample" = "sample","phenotype" = "phenotype", join_additional)

  samples_tidy %>% 
    inner_join(samples_tidy, by = join_by) %>%
    mutate(odds.x = (value.x / (1 - value.x)),
           odds.y = (value.y / (1 - value.y)), odds_ratio =  odds.x / odds.y)
}

get_tidy_samples_prediction <- function(fit, data_for_prediction) {
  fitted_detailed <- posterior_predict(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as.tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) 
  
}
