intercept_prior <- prior(normal(0,2), class = Intercept)
sd_prior <- prior(normal(0,2), class = sd)

phenotypes_to_use <- c("RD","OBE","PD","CI","REN","LIV","REPROD","HEART")

phenotypes_to_show <- phenotypes_to_use

genes_to_show_from_data <- function(data) {
  data %>% select(gene,functional_group) %>% distinct() %>% 
    filter(functional_group != c("Others")) %>% get("gene", .)
} 


formula_gene_only <- brmsformula(phenotype_value ~ (1||phenotype) + ((0 + phenotype)|gene) , family = "bernoulli")

formula_gene_source <- brmsformula(update(formula_gene_only, ~ . + ((0 + phenotype)||source)) , family = "bernoulli")
