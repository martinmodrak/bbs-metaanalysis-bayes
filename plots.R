posterior_interval <- 0.95
posterior_interval_bounds <- c(0.5 * (1 - posterior_interval), 1 - (0.5 * (1 - posterior_interval)))


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
                          prediction_filter = NULL,
                          out_func = print) {
  
  predicted <- posterior_predict(fit, nsamples = 1000)
  data <- data_long
  
  if(!is.null(prediction_filter)) {
    predicted <- predicted[,prediction_filter] 
    data <- data_long[prediction_filter,]
  }

  observed <- data$phenotype_value
  
  if("overall" %in% types) {
    my_ppc_bars(data_filtered$phenotype_value, predicted) %>% out_func
  }
  if("gene" %in% types) {
    my_ppc_bars_grouped(data_filtered$phenotype_value, predicted, group = data$gene) %>% out_func
  }
  if("phenotype" %in% types) {
    my_ppc_bars_grouped(data_filtered$phenotype_value, predicted, group = data$phenotype) %>% out_func
  }
  
  if("functional_group" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = data$functional_group) %>% out_func
  }
  
  if("functional_group_phenotype" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = interaction(data$functional_group, data$phenotype)) %>% out_func
  }
  
  sex_fct <- factor(if_else(is.na(data$Sex), "NA", as.character(data$Sex)))
  if("sex" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = sex_fct) %>% out_func
  }
  
  if("sex_phenotype" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = interaction(sex_fct, data$phenotype)) %>% out_func
  }
  
  
  age_fct <- factor(if_else(is.na(data$age_group), "NA", as.character(data$age_group)))
  if("age" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = age_fct) %>% out_func
  }
  
  if("age_phenotype" %in% types) {
    my_ppc_bars_grouped(observed, predicted, group = interaction(age_fct, data$phenotype)) %>% out_func
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
    facet_wrap(~phenotype, ncol = 4)  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
  
}


get_tidy_samples <- function(fit, data_for_prediction) {
  fitted_detailed <- fitted(fit, data_for_prediction, summary = FALSE, allow_new_levels = TRUE, nsamples = 1000)
  
  samples_tidy <- data_for_prediction %>% 
    cbind(as.tibble(t(fitted_detailed))) %>%
    gather("sample","value", V1:V1000) %>%
    mutate(odds = value/(1-value))
  
}


simple_num_format <- function(x) {format(x, digits = 1, scientific = FALSE, drop0trailing = TRUE)}


plot_gene_phenotype_differences_estimates <- function(fit, data_for_prediction, genes_to_show = unique(data_for_prediction$gene)) {
  if(is.null(data_for_prediction$group)) {
    data_for_prediction$group <- data_for_prediction$gene
  }
  samples_tidy <- get_tidy_samples(fit, data_for_prediction)
  
  per_phenotype_and_sample_average <- samples_tidy %>%
    group_by(sample, phenotype) %>% 
    summarise(avg = mean(value), avg_odds = avg/(1-avg))
  
  data_to_plot <- samples_tidy %>% 
    filter(gene %in% genes_to_show) %>%
    inner_join(per_phenotype_and_sample_average, by = c("phenotype" = "phenotype", "sample" = "sample")) %>%
    #mutate(relative = value - avg) %>%
    mutate(relative = odds / avg_odds) %>%
    group_by(phenotype, group, functional_group) %>%
    summarise(Estimate = median(relative), 
              lower = quantile(relative, posterior_interval_bounds[1]),
              upper = quantile(relative, posterior_interval_bounds[2]),
              lower50 = quantile(relative, 0.25),
              upper50 = quantile(relative, 0.75)
    )
  
  data_to_plot %>% ggplot(aes(x = group, y = Estimate, ymin = lower, ymax = upper, color = functional_group)) +
    geom_hline(yintercept = 1, color = "darkred")+ 
    geom_linerange(aes(ymin = lower50, ymax = upper50), size = 2) +
    geom_linerange() + 
    facet_wrap(~phenotype, scales = "free_y", ncol = 4)  + 
    scale_y_log10("Odds ratio", labels = simple_num_format) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    guides(color = FALSE)
}

plot_pairwise_differences <- function(fit, data_for_prediction, 
                                      plot_types = c("heatmap_min","heatmap_max","heatmap_ci_excl","linerange_all"), 
                                      out_func = function(name, plot) {print(plot)},
                                      title_add = "") {
  samples_tidy <- get_tidy_samples(fit, data_for_prediction) 
  samples_diff <- samples_tidy %>% 
    select(-source) %>%
    inner_join(samples_tidy %>% select(-source), by = c("sample" = "sample","phenotype" = "phenotype")) %>%
    mutate(odds_ratio = (value.x / (1 - value.x)) / (value.y / (1 - value.y)))
  
  
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
    summarise(Estimate = median(odds_ratio), 
              lower = quantile(odds_ratio, posterior_interval_bounds[1]),
              upper = quantile(odds_ratio, posterior_interval_bounds[2]),
              lower50 = quantile(odds_ratio, 0.25),
              upper50 = quantile(odds_ratio, 0.75),
              one_location = ecdf(odds_ratio)(1),
              p_positive = mean(odds_ratio > 1),
              p_negative = mean(odds_ratio < 1)) %>%
    ungroup() %>%
    filter(gene.x != gene.y) %>%
    mutate(min_probable_OR = case_when(
      
        ((lower > 1) != (upper > 1)) ~ 1,
        lower < 1 ~ upper,
        lower > 1 ~ lower,
        TRUE ~ NA_real_ #Should not happen
      ),
      max_probable_OR = pmax(lower, upper, 1/lower, 1/upper),
      CI_excl_one = abs(one_location - 0.5) * 2,
      CI_excl_one_sign = if_else(p_positive > p_negative, CI_excl_one, - CI_excl_one),
      p_diff = if_else(p_positive > p_negative, p_positive, -p_negative),
      result_category = case_when(
        ((lower > 1) == (upper > 1)) ~ "95_excl",
        ((lower50 > 1) == (upper50 > 1)) ~ "50_excl",
        TRUE ~ "none"
      ),
      result_category = factor(result_category, levels = c("none","50_excl", "95_excl"))
    )
  
  theme_heatmap <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5), legend.text = NULL)

  if("heatmap_p" %in% plot_types) {
    lims <- c(-1, 1) 
    p <- data_for_diff  %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = p_diff)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradientn(limits = lims, 
                           colours = c("#ca0020","#f4a582","#f7f7f7","#f7f7f7","#92c5de","#0571b0"),
                           values = (c(-1,-0.75,-0.5,0.5,0.75,1) + 1) / 2) +
      theme_heatmap +
      facet_wrap(~phenotype, ncol = 4) +
      ggtitle(paste0("Probability the difference is systematic", title_add))
    out_func(paste0("heatmap_p", title_add), p)
  }
  
  if("heatmap_ci_excl" %in% plot_types) {
    lims <- c(-1, 1) 
    p <- data_for_diff  %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = CI_excl_one_sign)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradient2(limits = lims, 
                           low = "#ca0020", mid = "#f7f7f7", high = "#0571b0") +
      theme_heatmap +
      facet_wrap(~phenotype, ncol = 4) +
      ggtitle(paste0("Widest credible interval excluding 0", title_add))
    out_func(paste0("heatmap_ci_excl", title_add), p)
  }
  
  
  if("heatmap_min" %in% plot_types) {

    p <- data_for_diff %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = min_probable_OR)) +
        geom_tile(width = 1, height = 1) +
        scale_fill_gradient2(low = "#ca0020", high = "#0571b0", mid = "#f7f7f7", 
                             trans = "log10", breaks = c(0.5,1,2), limits = c(0.5,2)) +
        theme_heatmap +
        facet_wrap(~phenotype, ncol = 4) +
      ggtitle(paste0("Odds ratio, 95% conservative", title_add))
    out_func(paste0("heatmap_min", title_add), p)
    
  }
  

  
  if("heatmap_max" %in% plot_types) {
    p <- data_for_diff %>%  
      ggplot(aes(x = gene.x, y = gene.y, fill = max_probable_OR)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_distiller(palette = "Spectral", trans = "log10") +
      theme_heatmap +
      facet_wrap(~phenotype, ncol = 4)+
      ggtitle(paste0("Odds ratio, 95% maximal", title_add))
    out_func(paste0("heatmap_max", title_add), p)
  }
  
  
  scale_linerange <- scale_x_log10(labels = simple_num_format)
  scale_linerange_color <- scale_color_gradientn(limits = c(0,1), colours = c("#303030","#303030","#0571b0","#ca0020","#ca0020"),
                                                 values = c(0,0.4,0.6,0.95, 1), breaks = c(0,0.5,0.95))

  if("linerange" %in% plot_types) {
    for(ph in unique(data_for_prediction$phenotype)) {
      p <- data_for_diff %>% filter(phenotype == ph) %>%
        ggplot(aes(color = CI_excl_one)) + 
        geom_vline(xintercept = 1, color = "darkred") +
        # geom_segment(aes(x = lower, xend = upper), y = 0.5, yend = 0.5) + 
        # geom_segment(aes(x = lower50, xend = upper50), y = 0.5, yend = 0.5, size = 2) +
        geom_segment(aes(x = lower, xend = upper, y = gene.x, yend = gene.x)) + 
        geom_segment(aes(x = lower50, xend = upper50, y = gene.x, yend = gene.x), size = 2) +
        scale_linerange +
        scale_linerange_color +
        facet_grid( ~ gene.y)  +
        ggtitle(paste0("95% and 50% credible intervals for pairwise differences",ph , title_add))
      out_func(paste0("linerange_", ph, title_add), p)
    }
  }

  theme_linerange_all <- theme(axis.text=element_text(size=8), strip.text = element_text(size = 9))
  if("linerange_all" %in% plot_types) {
      p <- data_for_diff %>% 
        ggplot(aes(color = CI_excl_one)) + 
        geom_vline(xintercept = 1, color = "darkred") +
        # geom_segment(aes(x = lower, xend = upper), y = 0.5, yend = 0.5) + 
        # geom_segment(aes(x = lower50, xend = upper50), y = 0.5, yend = 0.5, size = 2) +
        geom_segment(aes(x = lower, xend = upper, y = gene.y, yend = gene.y)) + 
        geom_segment(aes(x = lower50, xend = upper50, y = gene.y, yend = gene.y), size = 2) +
        scale_linerange +
        scale_linerange_color +
        facet_grid(phenotype ~ gene.x)  +
        theme_linerange_all +
        ggtitle(paste0("95% and 50% credible intervals for pairwise differences", title_add))
      out_func(paste0("linerange_all", title_add), p)
  }
  
    
  if("linerange_all2" %in% plot_types) {
    p <- data_for_diff %>% 
      ggplot(aes(y = phenotype, yend = phenotype, color = CI_excl_one)) + 
      geom_vline(xintercept = 0, color = "darkred") +
      geom_segment(aes(x = lower, xend = upper)) + 
      geom_segment(aes(x = lower50, xend = upper50), size = 2) +
      scale_linerange +
      scale_linerange_color +
      facet_grid(gene.x ~ gene.y) +
      theme_linerange_all  +
      ggtitle(paste0("95% and 50% credible intervals for pairwise differences", title_add))
    out_func(paste0("linerange_all2", title_add), p)
    out_func(p)
  }  
}