
phenotypes_to_use <- c("RD","OBE","PD","CI","REP","REN","HEART","LIV","DD")
phenotypes_to_use_factor <- factor(phenotypes_to_use, levels = phenotypes_to_use)

phenotypes_to_show <- phenotypes_to_use_factor

genes_to_show_from_data <- function(data) {
  data %>% select(gene,functional_group) %>% distinct() %>% 
    filter(functional_group != c("Others")) %>% get("gene", .)
}


intercept_prior <- prior(normal(0,2), class = Intercept)
b_prior <- prior(normal(0,2), class = b)
sd_prior <- prior(normal(0,2), class = sd)
family_sd_prior <- prior(normal(0,1), class = sd, group = "family_id")

default_priors <- c(intercept_prior, sd_prior)

family_default_priors <- c(default_priors, family_sd_prior)

formula_gene_only <- phenotype_value ~ (1||phenotype) + ((0 + phenotype)|gene) 

formula_gene_source <- update(formula_gene_only, ~ . + ((0 + phenotype)||source)) 

add_lof <- function(original_formula) {
  update(original_formula, ~ . +  phenotype : loss_of_function_certain) 
}
add_lof_per_gene <- function(original_formula) {
  update(original_formula, ~ . + ((0 + phenotype : loss_of_function_certain)||gene)) 
}

add_family <- function(orginal_formula) {
  update.formula(orginal_formula,  ~ . + ((0 + phenotype)||family_id))
}

add_ethnic_group <- function(original_formula) {
  update.formula(original_formula,  ~ . + ((0 + phenotype)||ethnic_group))
}

add_ethnicity <- function(original_formula) {
  update.formula(original_formula,  ~ . + ((0 + phenotype)||ethnicity))
}

add_sex <- function(original_formula) {
  update.formula(original_formula,  ~ . + ((0 + phenotype)||sex))
}

add_age <- function(original_formula) {
  update.formula(original_formula, ~ . + (0 + age_std_for_model||phenotype))
}

formula_gene_lof <- add_lof(formula_gene_only)
formula_gene_lof_per_gene <- add_lof_per_gene(formula_gene_lof)

formula_gene_family <- add_family(formula_gene_only)
formula_gene_ethnic_group <- add_ethnic_group(formula_gene_only)
formula_gene_ethnicity <- add_ethnicity(formula_gene_only)

formula_gene_ethnic_group_lof <- add_lof(formula_gene_ethnic_group)
formula_gene_ethnicity_lof <- add_lof(formula_gene_ethnicity)

formula_gene_ethnic_group_lof_per_gene <- add_lof_per_gene(formula_gene_ethnic_group_lof)
formula_gene_ethnicity_lof_per_gene <- add_lof(formula_gene_ethnicity_lof)

formula_gene_family_age_sex <- add_age(add_sex(formula_gene_family))

formula_gene_source_lof <- add_lof(formula_gene_source)
formula_gene_source_lof_per_gene <- add_lof_per_gene(formula_gene_source_lof)

formula_gene_source_family <- add_family(formula_gene_source)
formula_gene_source_lof_family <- add_family(formula_gene_source_lof)

formula_gene_source_lof_ethnic_group <- add_ethnic_group(formula_gene_source_lof)

formula_gene_source_sex <- add_sex(formula_gene_source)

formula_gene_source_age <- add_age(formula_gene_source)                                       

formula_gene_age_sex <- add_age(formula_gene_source_sex) 

formula_gene_age_sex_lof <- add_lof(formula_gene_age_sex)
formula_gene_age_sex_lof_per_gene <- add_lof_per_gene(formula_gene_age_sex_lof)

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
    formula = formula_gene_lof, 
    priors = c(default_priors, b_prior), 
    note = "",
    additional_components = c("lof"), filter = "none"
  ),
  gene_lof_per_gene = list(
    formula = formula_gene_lof_per_gene, 
    priors = c(default_priors, b_prior), 
    note = "",
    additional_components = c("lof per gene"), filter = "none"
  ),
  
  gene_family = list(
    formula = formula_gene_family, 
    priors = family_default_priors, 
    note = "",
    additional_components = c("family"), filter = "family"
  ),

  gene_ethnic_group = list(
    formula = formula_gene_ethnic_group, 
    priors = default_priors, 
    note = "",
    additional_components = c("ethnic_group"), filter = "ethnic_group"
  ),

  gene_ethnicity = list(
    formula = formula_gene_ethnicity, 
    priors = default_priors, 
    note = "",
    additional_components = c("ethnicity"), filter = "ethnicity"
  ),
  
  gene_ethnic_group_lof = list(
    formula = formula_gene_ethnic_group_lof, 
    priors = default_priors, 
    note = "",
    additional_components = c("lof", "ethnic_group"), filter = "ethnic_group"
  ),

  gene_ethnicity_lof = list(
    formula = formula_gene_ethnicity_lof, 
    priors = default_priors, 
    note = "",
    additional_components = c("lof","ethnicity"), filter = "none"
  ),
  
  gene_ethnic_group_lof_per_gene = list(
    formula = formula_gene_ethnic_group_lof_per_gene, 
    priors = default_priors, 
    note = "",
    additional_components = c("lof per gene","ethnic_group"), filter = "none"
  ),

  gene_ethnicity_lof_per_gene = list(
    formula = formula_gene_ethnicity_lof_per_gene, 
    priors = default_priors, 
    note = "",
    additional_components = c("lof per gene","ethnicity"), filter = "none"
  ),
  
  gene_family_filtered_age_sex_lof = list(
    formula = formula_gene_family_age_sex, 
    priors = family_default_priors, 
    note = "",
    additional_components = c("family", "age", "sex"), filter = c("family","age","sex","lof")
  ),
  
  gene_source = list(
    formula = formula_gene_source, 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "none"
    ),
  

  gene_source_lof =  list(
    formula = formula_gene_source_lof, 
    priors = c(default_priors, b_prior), 
    note = "",
    additional_components = c("source", "lof"), filter = "none"
  ),
  gene_source_family = list(
    formula = formula_gene_source_family, 
    priors = family_default_priors, 
    note = "",
    additional_components = c("source","family"), filter = "family"
  ),
  gene_source_lof_family =  list(
    formula = formula_gene_source_lof_family,
    priors = c(default_priors, b_prior, family_sd_prior),
    note = "",
    additional_components = c("source", "lof", "family"), filter = "family"
  ),
  gene_source_lof_ethnic_group =  list(
    formula = formula_gene_source_lof_ethnic_group,
    priors = c(default_priors, b_prior),
    note = "",
    additional_components = c("source", "lof", "ethnic_group"), filter = "ethnic_group"
  ),
  
  gene_source_lof_wide = list(
    formula = formula_gene_source_lof, 
    priors = c(prior(normal(0,5), class = Intercept), prior(normal(0,5), class = sd), b_prior), 
    note = "special priors",
    additional_components = c("source","lof"), filter = "none"
  ),
  

  gene_source_filtered_lof = list(
    formula = formula_gene_source, 
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "lof"
  ),

  gene_source_filtered_lof_wide = list(
    formula = formula_gene_source, 
    priors = c(prior(normal(0,5), class = Intercept), prior(normal(0,5), class = sd)), 
    note = "special priors",
    additional_components = c("source"), filter = "lof"
  ),

  gene_source_lof_per_gene =  list(
    formula = formula_gene_source_lof_per_gene, 
    priors = c(default_priors, b_prior), 
    note = "",
    additional_components = c("source", "lof per gene"), filter = "none"
  ),
  
  gene_source_genecor = list(
    formula = phenotype_value ~ (1 || phenotype) + ((0 + gene) | phenotype) + 
                            ((0 + phenotype) || source), 
                           
    priors = default_priors, 
    note = "",
    additional_components = c("source"), filter = "none"
  ),
  
  gene_source_nocor = list(
    formula = phenotype_value ~ (1 || phenotype) + ((0 + phenotype)||gene) + 
                            ((0 + phenotype) || source)  , 
                           
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
    additional_components = c("age","sex"), filter = c("age","sex")
  ),
  
  gene_lof_filtered_age_sex = list(
    formula = formula_gene_age_sex_lof, 
    priors = default_priors, 
    note = "",
    additional_components = c("age","sex", "lof"), filter = c("age","sex")
  ),
  gene_lof_per_gene_filtered_age_sex = list(
    formula = formula_gene_age_sex_lof_per_gene, 
    priors = c(default_priors, b_prior), 
    note = "",
    additional_components = c("age","sex", "lof per gene"), filter = c("age","sex")
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
    formula = update.formula(formula_gene_source, ~ . + (0 + age_std_for_model||phenotype)), 
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

main_model_def <- models_base$gene_source_lof
