
phenotypes_to_use <- c("RD","OBE","PD","CI","REP","REN","HEART","LIV")
phenotypes_to_use_factor <- factor(phenotypes_to_use, levels = phenotypes_to_use)

phenotypes_to_show <- phenotypes_to_use_factor

genes_to_show_from_data <- function(data) {
  data %>% select(gene,functional_group) %>% distinct() %>% 
    filter(functional_group != c("Others")) %>% get("gene", .)
} 


intercept_prior <- prior(normal(0,2), class = Intercept)
sd_prior <- prior(normal(0,2), class = sd)

default_priors <- c(intercept_prior, sd_prior)


formula_gene_only <- brmsformula(phenotype_value ~ (1||phenotype) + ((0 + phenotype)|gene) , family = "bernoulli")

formula_gene_source <- brmsformula(update(formula_gene_only, ~ . + ((0 + phenotype)||source)) , family = "bernoulli")

formula_gene_source_lof <- brmsformula(update(formula_gene_source, ~ . + ((0 + phenotype : loss_of_function_certain)||gene)) , family = "bernoulli")


formula_gene_source_sex <- brmsformula(update.formula(formula_gene_source,  ~ . + (1||Sex : phenotype)), family = "bernoulli")

formula_gene_source_age <- brmsformula(update.formula(formula_gene_source, ~ . + (0 + age_std_for_model||phenotype)), family = "bernoulli")                                       

formula_gene_age_sex <- brmsformula(update(formula_gene_only, ~ . + (1||Sex : phenotype) + (0 + age_std_for_model||phenotype)) , 
                               family = "bernoulli")

models_base <- list(
  gene_only = list(
    formula = formula_gene_only, 
    priors = default_priors, 
    note = "",
    additional_components = c(), filter = "none"
    ),
  
  gene_only_filtered_lof = list(
    formula = formula_gene_only, 
    priors = default_priors, 
    note = "",
    additional_components = c(), filter = "lof"
    ),
  
  gene_lof = list(
    formula = brmsformula(update(formula_gene_only, ~ . + ((0 + phenotype : loss_of_function_certain)||gene)) , 
                          family = "bernoulli"), 
    priors = default_priors, 
    note = "",
    additional_components = c("lof"), filter = "none"
  ),

  gene_source = list(
    formula = formula_gene_source, 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "none"
    ),
  
  gene_source_filtered_lof = list(
    formula = formula_gene_source, 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "lof"
    ),

  gene_source_lof =  list(
    formula = formula_gene_source_lof, 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "lof"), filter = "none"
  ),
  
  gene_source_genecor = list(
    formula = brmsformula(phenotype_value ~ (1 || phenotype) + ((0 + gene) | phenotype) + 
                            ((0 + phenotype) || source), 
                          family = "bernoulli"), 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "none"
  ),
  
  gene_source_nocor = list(
    formula = brmsformula(phenotype_value ~ (1 || phenotype) + ((0 + phenotype)||gene) + 
                            ((0 + phenotype) || source)  , 
                          family = "bernoulli"), 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "none"
  ),  
  
  gene_source_narrow = list(
    formula = formula_gene_source, 
    priors = c(intercept_prior, sd_prior, prior(normal(0, 1), class = sd, group = gene)), 
    note = "special priors",
    additional_components = c("source"), filter = "none"
  ),
  
  gene_source_very_narrow = list(
    formula = formula_gene_source, 
    priors = c(intercept_prior, sd_prior, prior(normal(0,0.1), class = sd, group = gene)), 
    note = "special priors",
    additional_components = c("source"), filter = "none"
  ),
  
  gene_source_wide = list(
    formula = formula_gene_source, 
    priors = c(prior(normal(0,5), class = Intercept), prior(normal(0,5), class = sd)), 
    note = "special priors",
    additional_components = c("source"), filter = "none"
  ),
  
  gene_source_filtered_sex = list(
    formula = formula_gene_source_sex, 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "sex"), filter = "sex"
    ),
  
  gene_source_filtered_age = list(
    formula = formula_gene_source_age, 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "age"), filter = "age"
    ),
  
  gene_filtered_age_sex = list(
    formula = formula_gene_age_sex, 
    priors = default_priors, 
    note = "",
    additional_components = c("age","sex"), filter = "age_sex"
  )
  
) %>% imap(
  function(def, name) { 
    def$name <- name
    def$imputed <- FALSE
    def
    })


models_imputed <- list(
  gene_imputed_age_sex = list(
    formula = formula_gene_age_sex, 
    priors = default_priors, 
    note = "",
    additional_components = c("age", "sex")
  ),
  
  gene_source_imputed_sex = list(
    formula = formula_gene_source_sex, 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "sex")
    ),
  
  gene_source_imputed_age = list(
    formula = formula_gene_source_age, 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "age")
    ),
  
  gene_source_imputed_age_sex = list(
    formula = brmsformula(update.formula(formula_gene_source, ~ . + (0 + age_std_for_model||phenotype)), 
                          family = "bernoulli"), 
    priors = default_priors, 
    note = "",
    additional_components = c("source", "age", "sex")
    )
) %>% imap(
  function(def, name) { 
    def$name <- name
    def$imputed <- TRUE
    def
  })

all_models <- c(models_base, models_imputed)