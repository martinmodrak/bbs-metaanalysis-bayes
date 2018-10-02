plot_gene_phenotype_estimates <- function(fit, data_for_prediction) {
  estimates <- fitted(fit, data_for_prediction, allow_new_levels = TRUE, probs = posterior_interval_bounds, robust = TRUE)
  estimates50 <- fitted(fit, data_for_prediction, allow_new_levels = TRUE, probs = c(0.25,0.75), robust = TRUE)
  
  
  data_with_prediction <- data_for_prediction %>%
    cbind(as.tibble(estimates)) %>% cbind(as.tibble(estimates50) %>% select(Q25,Q75)) 
  
  data_with_prediction$lower <- data_with_prediction[,paste0("Q",posterior_interval_bounds[1] * 100)]
  data_with_prediction$upper <- data_with_prediction[,paste0("Q",posterior_interval_bounds[2] * 100)]
  
  data_with_prediction %>% ggplot(aes(x = gene, y = Estimate, ymin = lower, ymax = upper, color = gene)) +
    geom_linerange() + 
    geom_linerange(aes(ymin = Q25, ymax = Q75), size = 2) +
    facet_wrap(~phenotype)  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
  
}


get_tidy_samples <- function(fit, data_for_prediction) {
  fitted_detailed <- fitted(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as.tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V4000) 
  
}


plot_gene_phenotype_differences_estimates <- function(fit, data_for_prediction, genes_to_show = unique(data_for_prediction$gene)) {
  samples_tidy <- get_tidy_samples(fit, data_for_prediction)
  
  per_phenotype_and_sample_average <- samples_tidy %>%
    group_by(sample, phenotype) %>% 
    summarise(avg = mean(value))
  
  data_to_plot <- samples_tidy %>% 
    filter(gene %in% genes_to_show) %>%
    inner_join(per_phenotype_and_sample_average, by = c("phenotype" = "phenotype", "sample" = "sample")) %>%
    mutate(relative = value - avg) %>%
    group_by(phenotype, gene) %>%
    summarise(Estimate = median(relative), 
              lower = quantile(relative, posterior_interval_bounds[1]),
              upper = quantile(relative, posterior_interval_bounds[2]),
              lower50 = quantile(relative, 0.25),
              upper50 = quantile(relative, 0.75)
    )
  
  data_to_plot %>% ggplot(aes(x = gene, y = Estimate, ymin = lower, ymax = upper, color = gene)) +
    geom_hline(yintercept = 0, color = "darkred")+ 
    geom_linerange(aes(ymin = lower50, ymax = upper50), size = 2) +
    geom_linerange() + facet_wrap(~phenotype, scales = "free_y")  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
}

plot_pairwise_differences <- function(fit, data_for_prediction) {
  
}