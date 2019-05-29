posterior_interval <- 0.95
posterior_interval_bounds <- c(0.5 * (1 - posterior_interval), 1 - (0.5 * (1 - posterior_interval)))


base_theme_font_size = (9 + 1/3) * 1.4
base_theme <- theme( 
                    axis.title = element_text(size = base_theme_font_size), 
                    axis.text = element_text(size = base_theme_font_size), 
                    legend.title = element_text(size = base_theme_font_size + 1),
                    legend.text = element_text(size = base_theme_font_size),
                    strip.text = element_text(size = base_theme_font_size),
                    plot.title = element_text(size = base_theme_font_size + 2))

scale_color_functional_group <- scale_color_manual("Functional group", values = c(BBS03 = "#f8766dff", BBSome = "#5f8dd3ff", Chaperonins = "#aa8800ff", Others = "#43aa8b"))

bbs_labeller <- function(x) {
  gsub("BBS0","BBS", as.character(x), fixed = TRUE)
}

bbs_labeller_facet <- labeller(gene = bbs_labeller)


my_ppc <- function(fun, ...) {
  #I unfortunatel need different settigns from within RStudio and for PDF render
  if(is.null(knitr::current_input())) {
    fatten = 1.5
    size = 1
  } else {
    fatten = 0.8
    size = 0.5
  }
  fun(..., prob = posterior_interval, fatten = fatten, size = size, freq = FALSE)
}

my_ppc_bars <- function(...) {
  my_ppc(ppc_bars, ...) + theme(axis.text = element_blank())
}

my_ppc_bars_grouped <- function(group, ...) {
  n_groups <- length(unique(group))
  if(ceiling(n_groups / 7) < (n_groups / 6)) {
    ncol = 7
  } else {
    ncol = 6
  }
  res <- my_ppc(ppc_bars_grouped, facet_args = list(ncol = ncol), group = group, ...) + 
    theme(axis.text = element_blank())
  
  suppressMessages(res  +  scale_x_continuous(breaks = c(0,1)))
}

run_pp_checks <- function(model_def, fit, data_long, 
                          types = c("overall","gene","phenotype","functional_group","functional_group_phenotype",
                                    "sex", "sex_phenotype","age","age_phenotype", "gene_lof", "phenotype_lof"), 
                          prediction_filter = NULL,
                          out_func = print) {
  
  predicted <- posterior_predict(fit, nsamples = 1000)
  data <- filter_data_by_model_def(model_def, data_long)
  
  if(!is.null(prediction_filter)) {
    predicted <- predicted[,prediction_filter] 
    data <- data_long[prediction_filter,]
  }

  observed <- data$phenotype_value
  
  gene_groups = if_else(data$functional_group != "Others" & data$gene != "BBS18", as.character(data$gene), "other")
  
  small_labels <- theme(strip.text = element_text(size = 6))
  
  
  if("overall" %in% types) {
    (my_ppc_bars(data$phenotype_value, predicted) + ggtitle(paste0("PPCheck overall for ", model_def$name))) %>% out_func
  }
  if("gene" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = data$gene) + ggtitle(paste0("PPCheck gene for ", model_def$name))) %>% out_func
  }
  if("source" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = data$source) + 
       ggtitle(paste0("PPCheck source for ", model_def$name)) +
       small_labels
     ) %>% out_func
  }
  
  filter_10 <- data %>% group_by(source) %>% mutate(include = length(unique(ID)) >= 10) %>% pull(include)
  if("source_10" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value[filter_10], predicted[, filter_10], group = data$source[filter_10]) + 
       ggtitle(paste0("PPCheck sources with >= 10 patients for ", model_def$name)) +
       small_labels
    ) %>% out_func
  }
  if("source_10_lof" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value[filter_10], predicted[, filter_10], group = interaction(data$loss_of_function[filter_10], data$source[filter_10])) + 
       ggtitle(paste0("PPCheck sources with >= 10 patients for ", model_def$name)) +
       small_labels
    ) %>% out_func
  }
  
  filter_family_4 <- data %>% group_by(family_id) %>% mutate(include = !any(is.na(family_id)) & length(unique(ID)) >= 4) %>% pull(include)
  if("family_4" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value[filter_family_4], predicted[, filter_family_4], group = data$family_id[filter_family_4]) + 
       ggtitle(paste0("PPCheck families with >= 4 patients for ", model_def$name)) +
       small_labels
    ) %>% out_func
  }

  filter_family_3 <- data %>% group_by(family_id) %>% mutate(include = !any(is.na(family_id)) & length(unique(ID)) == 3) %>% pull(include)
  if("family_3" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value[filter_family_3], predicted[, filter_family_3], group = data$family_id[filter_family_3]) + 
       ggtitle(paste0("PPCheck families with exactly 3 patients for ", model_def$name)) +
       small_labels
    ) %>% out_func
  }
  if("ethnic_group" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = data$ethnic_group) + ggtitle(paste0("PPCheck ethnic group for ", model_def$name))) %>% out_func
  }
  
  filter_ethnicity_10 <- data %>% group_by(ethnicity) %>% mutate(include = length(unique(ID)) >= 10) %>% pull(include)
  
  if("ethnicity_10" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value[filter_ethnicity_10], predicted[, filter_ethnicity_10], group = data$ethnicity[filter_ethnicity_10]) + ggtitle(paste0("PPCheck ethnicities with >= 10 patients for ", model_def$name))) %>% out_func
  }
  
  
  if("gene_lof" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = interaction(data$loss_of_function, gene_groups)) + ggtitle(paste0("PPCheck gene_lof for ", model_def$name))) %>% out_func
  }
  if("phenotype" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = data$phenotype) + ggtitle(paste0("PPCheck phenotype for ", model_def$name))) %>% out_func
  }
  if("phenotype_lof" %in% types) {
    (my_ppc_bars_grouped(data$phenotype_value, predicted, group = interaction(data$loss_of_function, data$phenotype)) + ggtitle(paste0("PPCheck phenotype + LOF for ", model_def$name))) %>% out_func
  }
  
  
  if("functional_group" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = data$functional_group) + ggtitle(paste0("PPCheck functional group for ", model_def$name))) %>% out_func
  }
  
  if("functional_group_phenotype" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = interaction(data$functional_group, data$phenotype)) + ggtitle(paste0("PPCheck group + phenotype for ", model_def$name))) %>% out_func
  }
  
  sex_fct <- factor(if_else(is.na(data$sex), "NA", as.character(data$sex)))
  if("sex" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = sex_fct) + ggtitle(paste0("PPCheck sex for ", model_def$name))) %>% out_func
  }
  
  if("sex_phenotype" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = interaction(sex_fct, data$phenotype)) + ggtitle(paste0("PPCheck sex + phenotype for ", model_def$name))) %>% out_func
  }
  
  
  age_fct <- factor(if_else(is.na(data$age_group), "NA", as.character(data$age_group)))
  if("age" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = age_fct) + ggtitle(paste0("PPCheck age for ", model_def$name))) %>% out_func
  }
  
  if("age_phenotype" %in% types) {
    age_fct_collapsed <- fct_collapse(age_fct, `0-19` = c("0-9","10-19"), `20-39` = c("20-29","30-39"), `40+` = c("40-49","50-59","60+"))
    (my_ppc_bars_grouped(observed, predicted, group = interaction(age_fct_collapsed, data$phenotype)) + ggtitle(paste0("PPCheck age + phenotype for ", model_def$name))) %>% out_func
  }

  if("gene_phenotype" %in% types) {
    (my_ppc_bars_grouped(observed, predicted, group = interaction(gene_groups, data$phenotype)) + ggtitle(paste0("PPCheck gene + phenotype for ", model_def$name))) %>% out_func
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
    facet_wrap(~phenotype, ncol = 3)  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    base_theme
  
}


simple_num_format <- function(x) {format(x, digits = 1, scientific = FALSE, drop0trailing = TRUE)}


plot_gene_phenotype_differences_estimates <- function(
  fit, data_for_prediction, 
  genes_to_show = unique(data_for_prediction$gene), 
  group_title = "Gene", data_original = NULL,
  wider = FALSE) {
  if(is.null(data_for_prediction[["group"]])) {
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
    ) %>%
    mutate(phenotype_long = phenotype_long_from_phenotype(phenotype))

  

  if(is.null(data_original)) {
    original_geom = NULL 
    limits_geom = NULL
  } else {
    if(is.null(data_original[["group"]])) {
      data_original$group <- data_original$gene
    }
    
    data_original_matched <- data_original %>%
      filter(phenotype %in% unique(data_for_prediction$phenotype),
             gene %in% genes_to_show,
             gene %in% unique(data_for_prediction$gene))
    
    per_phenotype_and_source_average_original <- data_original_matched %>%
      group_by(phenotype, source, gene) %>% #Two step averaging to reflect that all genes have the same weight in fitted data
      summarise(avg_per_gene = mean(phenotype_value)) %>%
      group_by(phenotype, source) %>%
      summarise(avg = mean(avg_per_gene)) %>%
      mutate(avg_odds = avg / (1-avg))
    
    data_to_plot_original <- data_original_matched %>%
      group_by(phenotype, group, source) %>%
      summarise(avg = mean(phenotype_value), num_cases = length(phenotype_value)) %>%
      mutate(odds = avg / (1-avg)) %>%
      inner_join(per_phenotype_and_source_average_original, by = c("phenotype" = "phenotype", "source" = "source")) %>%
      mutate(odds_ratio = case_when(
        avg.y == 0 ~ 1,
        avg.y == 1 ~ 1,
        avg.x == 0 ~ NA_real_,
        avg.x == 1 ~ NA_real_,
        TRUE ~ odds / avg_odds) )
    
    phenotype_limits_points <- data_to_plot_original %>% 
      filter(!is.na(odds_ratio), odds_ratio > 0) %>%
      group_by(phenotype) %>%
      summarise(min_or_points = min(odds_ratio), max_or_points = max(odds_ratio))

    phenotype_limits_intervals <- data_to_plot %>% 
      group_by(phenotype) %>%
      summarise(min_or_intervals = min(lower), max_or_intervals = max(upper))
    
    phenotype_limits <- phenotype_limits_points %>%
      inner_join(phenotype_limits_intervals, by = c("phenotype" = "phenotype")) %>%
      mutate(min_or = pmin(min_or_intervals, min_or_points), 
             max_or = pmax(max_or_intervals, max_or_points),
             limits_ratio = max_or / min_or,
             phenotype_long = phenotype_long_from_phenotype(phenotype)
             )
    
    data_to_plot_original_imputed <- data_to_plot_original %>%
      inner_join(phenotype_limits, by = c("phenotype"= "phenotype")) %>%
      mutate(type = case_when(
        !is.na(odds_ratio) ~ "normal",
        avg.x == 0 ~ "extreme",
        avg.x == 1 ~ "extreme",
        TRUE ~ NA_character_
      ),
      odds_ratio = case_when(
        !is.na(odds_ratio) ~ odds_ratio,
        avg.x == 0 ~ min_or / (limits_ratio ^ 0.1),
        avg.x == 1 ~ max_or * (limits_ratio ^ 0.1),
        TRUE ~ NA_real_)
      )
    
    if(wider) {
      jitter_width = 0.3
    } else {
      jitter_width = 0.2
    }
    original_geom = geom_point(data = data_to_plot_original_imputed, 
                               aes(size = num_cases, x = group, y = odds_ratio), 
                               inherit.aes = FALSE, color = "darkgray", alpha = 0.5,
                               position = position_jitter(width = jitter_width, height = 0))
    
    limits_geom = geom_hline(
      data = phenotype_limits %>%
        mutate(min_or_bound = min_or / (limits_ratio ^ 0.05),
               max_or_bound = max_or * (limits_ratio ^ 0.05)) %>%
        gather("type","value", min_or_bound, max_or_bound),
      aes(yintercept = value), 
      color = "blue",
      linetype = "dashed"
      )
  }
  
  if(wider) {
    size_50 = 3
    size_95 = 1.5
    size_range = c(1,4)
    facet <- facet_wrap(~phenotype_long, scales = "free_y", ncol = 3)
  } else {
    size_50 = 2
    size_95 = 1
    size_range = c(0.5,3)
    facet <- facet_wrap(~phenotype, scales = "free_y", ncol = 3)
  }
    
  data_to_plot %>% ggplot(aes(x = group, y = Estimate, ymin = lower, ymax = upper, color = functional_group)) +
    geom_hline(yintercept = 1, color = "darkred")+ 
    limits_geom +
    original_geom +
    geom_linerange(aes(ymin = lower50, ymax = upper50), size = size_50) +
    geom_linerange() +
    facet +
    scale_color_functional_group +
    scale_y_log10("Odds ratio", labels = simple_num_format) +
    scale_x_discrete(group_title, labels = bbs_labeller) +
    scale_alpha_continuous(range = c(0.1,0.6)) +
    scale_size_continuous(range = size_range) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    base_theme +
    guides(color = FALSE, size = FALSE, alpha = FALSE)
}

plot_correlations <- function(model_def, fit, data_for_prediction, 
                              out_func = function(name, plot) {print(plot)},
                              title_add = "") {
  
  samples_tidy <- get_tidy_samples(fit, data_for_prediction) 
  samples_diff <- get_samples_pairwise_diff(model_def, samples_tidy)
  
  if("cor" == plot_type) {
    samples_to_show = sample(unique(samples_tidy$sample), 100)
    data_for_cor <- samples_diff %>%
      filter(sample %in% samples_to_show) 
    
    for(ph in unique(data_for_prediction$phenotype)) {
      p <- data_for_cor %>% filter(phenotype == ph) %>%
        ggplot(aes(x = odds.x, y = odds.y)) + 
        geom_point(size = 0.1) + facet_grid(gene.x ~ gene.y, scales = "free") +
        scale_x_log10() + scale_y_log10() +
        ggtitle(paste0("Correlations for ", ph)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
      out_func(paste0("cor_",ph), p)
    }
  }  
}

plot_pairwise_differences <- function(model_def, fit, data_for_prediction, 
                                      plot_type = "linerange_all", 
                                      title_add = "") {
  samples_tidy <- get_tidy_samples(fit, data_for_prediction) 
  samples_diff <- get_samples_pairwise_diff(model_def, samples_tidy)
  
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
  
  theme_heatmap <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5),
                         axis.title = element_blank())

  if("heatmap_p" == plot_type) {
    lims <- c(-1, 1) 
    p <- data_for_diff  %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = p_diff)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradientn(limits = lims, 
                           colours = c("#ca0020","#f4a582","#f7f7f7","#f7f7f7","#92c5de","#0571b0"),
                           values = (c(-1,-0.75,-0.5,0.5,0.75,1) + 1) / 2) +
      theme_heatmap +
      facet_wrap(~phenotype, ncol = 3) +
      ggtitle(paste0("Probability the difference is systematic", title_add))
  }
  
  if("heatmap_ci_excl" == plot_type) {
    lims <- c(-1, 1) 
    p <- data_for_diff  %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = CI_excl_one_sign)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradient2(limits = lims, 
                           low = "#ca0020", mid = "#f7f7f7", high = "#0571b0") +
      base_theme +
      theme_heatmap +
      scale_x_discrete(labels = bbs_labeller) +
      scale_y_discrete(labels = bbs_labeller) +
      facet_wrap(~phenotype, ncol = 3) +
      ggtitle(paste0("Widest credible interval excluding 0", title_add))
  }
  
  
  if("heatmap_min" == plot_type) {

    p <- data_for_diff %>% 
      ggplot(aes(x = gene.x, y = gene.y, fill = min_probable_OR)) +
        geom_tile(width = 1, height = 1) +
        scale_fill_gradient2("     ", low = "#ca0020", high = "#0571b0", mid = "#f7f7f7", 
                             trans = "log10", breaks = c(0.5,1,2)) +
        expand_limits(fill = c( c(0.5,1,2))) +
        base_theme +
        theme_heatmap +
        scale_x_discrete(labels = bbs_labeller) +
        scale_y_discrete(labels = bbs_labeller) +
        facet_wrap(~phenotype, ncol = 3) +
      ggtitle(paste0("Odds ratio, 95% conservative", title_add))
  }
  

  
  if("heatmap_max" == plot_type) {
    p <- data_for_diff %>%  
      ggplot(aes(x = gene.x, y = gene.y, fill = max_probable_OR)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_distiller("     ", palette = "Spectral", trans = "log10") +
      base_theme +
      theme_heatmap +
      scale_x_discrete(labels = bbs_labeller) +
      scale_y_discrete(labels = bbs_labeller) +
      facet_wrap(~phenotype, ncol = 3)+
      ggtitle(paste0("Odds ratio, 95% extreme", title_add))
  }
  
  
  scale_linerange <- scale_x_log10(labels = simple_num_format)
  scale_linerange_color <- scale_color_gradientn("CI excluding 1" , limits = c(0,1), colours = c("#303030","#303030","#0571b0","#ca0020","#ca0020"),
                                                 values = c(0,0.4,0.6,0.95, 1), breaks = c(0,0.5,0.95),
                                                 labels = c("0%","50%","95%"))

  if("linerange" == plot_type) {
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
        base_theme +
        facet_grid( ~ gene.y)  +
        ggtitle(paste0("95% and 50% credible intervals for pairwise odds ratios",ph , title_add))
    }
  }

  # theme_linerange_all <- theme(axis.text=element_text(size=8), strip.text = element_text(size = 9),
  #                              axis.title = element_blank(), 
  #                              legend.title = element_text(size = 11,angle = 270, vjust = 0.2),
  #                              legend.text = element_text(size = 8))
  theme_linerange_all <- theme(axis.title = element_blank(), 
                               legend.title = element_text(angle = 270, vjust = 0.2))
  
  if("linerange_all" == plot_type) {
      p <- data_for_diff %>% 
        ggplot(aes(color = CI_excl_one)) + 
        geom_vline(xintercept = 1, color = "darkred") +
        # geom_segment(aes(x = lower, xend = upper), y = 0.5, yend = 0.5) + 
        # geom_segment(aes(x = lower50, xend = upper50), y = 0.5, yend = 0.5, size = 2) +
        geom_segment(aes(x = lower, xend = upper, y = gene.y, yend = gene.y)) + 
        geom_segment(aes(x = lower50, xend = upper50, y = gene.y, yend = gene.y), size = 2) +
        scale_linerange +
        scale_linerange_color +
        scale_y_discrete(labels = bbs_labeller) +
        facet_grid(phenotype ~ gene.x, labeller = labeller(gene.x = bbs_labeller))  +
        base_theme +
        theme_linerange_all +
        ggtitle(paste0("95% and 50% credible intervals for pairwise odds ratios", title_add))
  }
  
    
  if("linerange_all2" == plot_type) {
    p <- data_for_diff %>% 
      ggplot(aes(y = phenotype, yend = phenotype, color = CI_excl_one)) + 
      geom_vline(xintercept = 0, color = "darkred") +
      geom_segment(aes(x = lower, xend = upper)) + 
      geom_segment(aes(x = lower50, xend = upper50), size = 2) +
      scale_linerange +
      scale_linerange_color +
      facet_grid(gene.x ~ gene.y, labeller = labeller(gene.x = bbs_labeller, gene.y = bbs_labeller)) +
      theme_linerange_all  +
      base_theme +
      ggtitle(paste0("95% and 50% credible intervals for pairwise differences", title_add))
  }
  
  p
}


plot_lof_differences_estimates <- function(fit, data_for_prediction, genes_to_show = unique(data_for_prediction$gene), group_title = "Gene") {
  
  data_to_plot <- get_samples_lof_diff(fit, data_for_prediction) %>%
    group_by(phenotype, group, functional_group) %>%
    summarise(Estimate = median(odds_ratio), 
              lower = quantile(odds_ratio, posterior_interval_bounds[1]),
              upper = quantile(odds_ratio, posterior_interval_bounds[2]),
              lower50 = quantile(odds_ratio, 0.25),
              upper50 = quantile(odds_ratio, 0.75)
    )
  
  data_to_plot %>% ggplot(aes(x = group, y = Estimate, ymin = lower, ymax = upper, color = functional_group)) +
    geom_hline(yintercept = 1, color = "darkred")+ 
    geom_linerange(aes(ymin = lower50, ymax = upper50), size = 2) +
    geom_linerange() +
    facet_wrap(~phenotype, scales = "free_y", ncol = 3)  + 
    scale_color_functional_group +
    scale_y_log10("Odds ratio", labels = simple_num_format) +
    scale_x_discrete(group_title, labels = bbs_labeller) +
    scale_alpha_continuous(range = c(0.1,0.6)) +
    scale_size_continuous(range = c(0.5,3)) +
    base_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) +
    guides(color = FALSE, size = FALSE, alpha = FALSE)   +
    ggtitle("Odds ratio certain loss of function vs. others")
}
  