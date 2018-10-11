my_ppc <- function(fun, ...) {
  fun(..., prob = posterior_interval, fatten = 1.5, freq = FALSE)
}

my_ppc_bars <- function(...) {
  my_ppc(ppc_bars, ...)
}

my_ppc_bars_grouped <- function(...) {
  my_ppc(ppc_bars_grouped, ...)
}

run_pp_checks <- function(fit, data_long, 
                          types = c("overall","gene","phenotype","functional_group","functional_group_phenotype",
                                    "sex", "sex_phenotype","age","age_phenotype"), 
                          out_func = print) {
  
  predicted <- posterior_predict(fit, n_samples = 1000) 
  
  
  if("overall" %in% types) {
    my_ppc_bars(data_long$phenotype_value, predicted) %>% out_func
  }
  if("gene" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = data_long$gene) %>% out_func
  }
  if("phenotype" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = data_long$phenotype) %>% out_func
  }
  
  if("functional_group" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = data_long$functional_group) %>% out_func
  }
  
  if("functional_group_phenotype" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = interaction(data_long$functional_group, data_long$phenotype)) %>% out_func
  }
  
  sex_fct <- factor(if_else(is.na(data_long$Sex), "NA", as.character(data_long$Sex)))
  if("sex" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = sex_fct) %>% out_func
  }
  
  if("sex_phenotype" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = interaction(sex_fct, data_long$phenotype)) %>% out_func
  }
  
  
  age_fct <- factor(if_else(is.na(data_long$age_group), "NA", as.character(data_long$age_group)))
  if("age" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = age_fct) %>% out_func
  }
  
  if("age_phenotype" %in% types) {
    my_ppc_bars_grouped(data_long$phenotype_value, predicted, group = interaction(age_fct, data_long$phenotype)) %>% out_func
  }
  
}


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
  fitted_detailed <- fitted(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as.tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) %>%
    mutate(odds = value/(1-value))
  
}


plot_gene_phenotype_differences_estimates <- function(fit, data_for_prediction, genes_to_show = unique(data_for_prediction$gene)) {
  samples_tidy <- get_tidy_samples(fit, data_for_prediction)
  
  per_phenotype_and_sample_average <- samples_tidy %>%
    group_by(sample, phenotype) %>% 
    summarise(avg = mean(value), avg_odds = avg/(1-avg))
  
  data_to_plot <- samples_tidy %>% 
    filter(gene %in% genes_to_show) %>%
    inner_join(per_phenotype_and_sample_average, by = c("phenotype" = "phenotype", "sample" = "sample")) %>%
    #mutate(relative = value - avg) %>%
    mutate(relative = odds / avg_odds) %>%
    group_by(phenotype, gene) %>%
    summarise(Estimate = median(relative), 
              lower = quantile(relative, posterior_interval_bounds[1]),
              upper = quantile(relative, posterior_interval_bounds[2]),
              lower50 = quantile(relative, 0.25),
              upper50 = quantile(relative, 0.75)
    )
  
  data_to_plot %>% ggplot(aes(x = gene, y = Estimate, ymin = lower, ymax = upper, color = gene)) +
    geom_hline(yintercept = 1, color = "darkred")+ 
    geom_linerange(aes(ymin = lower50, ymax = upper50), size = 2) +
    geom_linerange() + facet_wrap(~phenotype, scales = "free_y")  + 
    scale_y_log10("Odds ratio") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
}

plot_pairwise_differences <- function(fit, data_for_prediction, plot_types = c("cor","heatmap_min","heatmap_max","heatmap_p", "linerange","linerange_all","linerange_all2"), out_func = print) {
  samples_tidy <- get_tidy_samples(fit, data_for_prediction) 
  samples_diff <- samples_tidy %>% 
    select(-source) %>%
    inner_join(samples_tidy %>% select(-source), by = c("sample" = "sample","phenotype" = "phenotype")) %>%
    mutate(difference = value.x - value.y)
  
  
  if("cor" %in% plot_types) {
    samples_to_show = sample(unique(samples_tidy$sample), 100)
    data_for_cor <- samples_diff %>%
      filter(sample %in% samples_to_show) 
    
    for(ph in unique(data_for_prediction$phenotype)) {
      p <- data_for_cor %>% filter(phenotype == ph) %>%
          ggplot(aes(x = value.x, y = value.y)) + 
          geom_point(size = 0.1) + facet_grid(gene.x ~ gene.y, scales = "free")
      out_func(p)
    }
   }
  
  data_for_diff <- samples_diff %>% 
    group_by(gene.x, gene.y, phenotype) %>%
    summarise(Estimate = median(difference), 
              lower = quantile(difference, posterior_interval_bounds[1]),
              upper = quantile(difference, posterior_interval_bounds[2]),
              lower50 = quantile(difference, 0.25),
              upper50 = quantile(difference, 0.75),
              p_positive = mean(difference > 0),
              p_negative = mean(difference < 0)) %>%
    mutate(min_probable_difference = case_when(
        sign(lower) != sign(upper) ~ 0,
        lower < 0 ~ upper,
        lower >= 0 ~ lower,
        TRUE ~ NA_real_ #Should not happen
      ),
      max_probable_difference = pmax(abs(lower), abs(upper)),
      p_diff = if_else(p_positive > p_negative, p_positive, -p_negative)
    )

  if("heatmap_p" %in% plot_types) {
    lims <- c(-1, 1) 
    p <- data_for_diff  %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = p_diff)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradientn(limits = lims, 
                           colours = c("#ca0020","#f4a582","#f7f7f7","#f7f7f7","#92c5de","#0571b0"),
                           values = (c(-1,-0.75,-0.5,0.5,0.75,1) + 1) / 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
      facet_wrap(~phenotype)
    out_func(p)
  }
  
  if("heatmap_min" %in% plot_types) {

    p <- data_for_diff %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = min_probable_difference)) +
        geom_tile(width = 1, height = 1) +
        scale_fill_gradient2(low = "#ca0020", high = "#0571b0", mid = "#f7f7f7") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
        facet_wrap(~phenotype)
    out_func(p)
    
  }
  

  
  if("heatmap_max" %in% plot_types) {
    lims <- c(0, 1) #max(data_for_diff$max_probable_difference))
    p <- data_for_diff %>%  
      ggplot(aes(x = gene.x, y = gene.y, fill = max_probable_difference)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_distiller(palette = "Spectral", limits = lims) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
      facet_wrap(~phenotype)
    out_func(p)
  }
  
  if("linerange" %in% plot_types) {
    for(ph in unique(data_for_prediction$phenotype)) {
      p <- data_for_diff %>% filter(phenotype == ph) %>%
        ggplot() + 
        geom_vline(xintercept = 0, color = "darkred") +
        # geom_segment(aes(x = lower, xend = upper), y = 0.5, yend = 0.5) + 
        # geom_segment(aes(x = lower50, xend = upper50), y = 0.5, yend = 0.5, size = 2) +
        geom_segment(aes(x = lower, xend = upper, y = gene.x, yend = gene.x)) + 
        geom_segment(aes(x = lower50, xend = upper50, y = gene.x, yend = gene.x), size = 2) +
        facet_grid( ~ gene.y)  +
        ggtitle(ph)
      out_func(p)
    }
  }

  if("linerange_all" %in% plot_types) {
      p <- data_for_diff %>% 
        ggplot() + 
        geom_vline(xintercept = 0, color = "darkred") +
        # geom_segment(aes(x = lower, xend = upper), y = 0.5, yend = 0.5) + 
        # geom_segment(aes(x = lower50, xend = upper50), y = 0.5, yend = 0.5, size = 2) +
        geom_segment(aes(x = lower, xend = upper, y = gene.x, yend = gene.x)) + 
        geom_segment(aes(x = lower50, xend = upper50, y = gene.x, yend = gene.x), size = 2) +
        facet_grid(phenotype ~ gene.y)  
        
      out_func(p)
  
  }
  
    
  if("linerange_all2" %in% plot_types) {
    p <- data_for_diff %>% 
      ggplot(aes(y = phenotype, yend = phenotype, color = phenotype)) + 
      geom_vline(xintercept = 0, color = "darkred") +
      geom_segment(aes(x = lower, xend = upper)) + 
      geom_segment(aes(x = lower50, xend = upper50), size = 2) +
      facet_grid(gene.x ~ gene.y)
    out_func(p)
  }  
}