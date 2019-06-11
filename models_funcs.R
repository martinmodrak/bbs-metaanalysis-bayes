filter_data_by_model_def <- function(def, data_long) {
  ret <- data_long
  for(filter in def$filter){
    ret <- switch (filter,
            "none" = ret,
            "lof" = ret %>% filter(loss_of_function == "certain"),
            "family" = ret %>% filter(!is.na(family_id)),
            "age" = ret %>% filter(!is.na(age_std_for_model)),
            "sex" = ret %>% filter(!is.na(sex)),
            "ethnic_group" = ret %>% filter(ethnic_group != "NA"),
            "ethnicity" = ret %>% filter(ethnicity != "NA"),
            stop("Unrecognized filter")
    )  
  }
  ret
}


stored_fits_dir <- "stored_fits"

fit_base_model <- function(def, data_long, ...) {
  if(is.null(def)) {
    stop("Invalid def")
  }
  
  if(def$imputed) {
    stop("Cannot be used on models with imputed data")
  }
  
  if(!dir.exists(here(stored_fits_dir))) {
    dir.create(here(stored_fits_dir))
  }

  data <- filter_data_by_model_def(def, data_long)
  
  brm(formula = brmsformula(def$formula, family = "bernoulli"), data = data, prior = def$priors, file = here(stored_fits_dir,def$name),
      control = list(adapt_delta = 0.95), ...)   
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
  
  brm_multiple(formula = brmsformula(def$formula, family = "bernoulli"), data = data_imputed, prior = def$priors, file = here(stored_fits_dir,def$name),
      control = list(adapt_delta = 0.95))   
}

data_for_prediction_base_model <- function(def, genes_to_show, phenotypes_to_show, 
                                           age_transform = NULL, 
                                           ages = c(10,30,50), sexes = c("F","M"),
                                           ethnic_groups = "EG-D", 
                                           ethnicities = c("Saudi Arabian", "French", "Pakistani", "Spanish", "Turkish"),
                                           loss_of_function_certain = 1) {
  result <-  tibble(gene = genes_to_show) %>%
    crossing(tibble(phenotype = phenotypes_to_show)) %>% 
    mutate(functional_group = functional_group_for_gene(gene))
  
  if("source" %in% def$additional_components) {
    result$source = "new_source" #A new factor level
  }

  if("family" %in% def$additional_components) {
    result$family_id = "new_family" #A new factor level
  }
  
  if("lof" %in% def$additional_components || "lof per gene" %in% def$additional_components) {
    result <- result %>% crossing(tibble(loss_of_function_certain = loss_of_function_certain)) %>%
      mutate(loss_of_function = if_else(loss_of_function_certain == 1, "certain", "unknown"))
  }
  
  if("age" %in% def$additional_components) {
    if(is.null(age_transform)) {
      stop("When working with age, age_transform has to be given")
    }
    result <- result %>% crossing(tibble(age_std_for_model = age_transform(ages)))
  }
  
  if("sex" %in% def$additional_components) {
    result <- result %>% crossing(tibble(sex = sexes))
  }
  
  if("ethnic_group" %in% def$additional_components) {
    result <- result %>% crossing(tibble(ethnic_group = ethnic_groups))
  }

  if("ethnicity" %in% def$additional_components) {
    result <- result %>% crossing(tibble(ethnicity = ethnicities))
  }
  
  result
}

filter_for_BBSome <- function(data_for_prediction) {
  filter(data_for_prediction, functional_group == "BBSome", gene != "BBS18")
}

get_tidy_samples <- function(fit, data_for_prediction, scale = "response") {
  fitted_detailed <- fitted(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000, scale = scale)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as_tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) 
  
  
  if(scale == "response") {
    samples_tidy %>%
      mutate(odds = value/(1-value))
  } else {
    samples_tidy
  }
  
}

get_samples_pairwise_diff <- function(model_def, samples_tidy) {
  component_to_join_map <- function(component) {
    switch(component,
           "lof" = c("loss_of_function_certain" = "loss_of_function_certain"),
           "age" = c("age_std_for_model" = "age_std_for_model"),
           "sex" = c("sex" = "sex"),
           "ethnic_group" = c("ethnic_group" = "ethnic_group"),
           "ethnicity" = c("ethnicity" = "ethnicity"),
           "family" = c("family_id" = "family_id"), 
           "source" = c("source" = "source")
    )
  }
  
  join_additional <- model_def$additional_components %>% map(component_to_join_map) %>% do.call(c, .)
  join_by <- c("sample" = "sample","phenotype" = "phenotype", join_additional)

  samples_tidy %>% 
    inner_join(samples_tidy, by = join_by) %>%
    mutate(odds_ratio =  odds.x / odds.y)
}

get_tidy_samples_prediction <- function(fit, data_for_prediction) {
  fitted_detailed <- posterior_predict(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as_tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) 
  
}

get_samples_lof_diff <- function(fit, data_for_prediction) {
  if(is.null(data_for_prediction[["group"]])) {
    data_for_prediction$group <- data_for_prediction$gene
  }
  
  data_for_prediction_lof <-  data_for_prediction %>% mutate(loss_of_function_certain = 1) %>% rbind(
    data_for_prediction %>% mutate(loss_of_function_certain = 0)
  )
  
  samples_tidy_lof <- get_tidy_samples(fit, data_for_prediction_lof)
  
  
  samples_tidy_lof %>% 
    filter(loss_of_function_certain == 1) %>%
    inner_join(samples_tidy_lof %>% filter(loss_of_function_certain == 0), 
               by = c("phenotype" = "phenotype", "gene" = "gene", "group" = "group", "sample" = "sample", 
                      "functional_group" = "functional_group")) %>%
    mutate(odds_ratio = odds.x / odds.y)  
}
