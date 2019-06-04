bbsome_genes <- paste0("BBS0", c(1,2,4,5,7, 8,9))

convert_uncertain_phenotype <- function(column) {
  result <- as.numeric(gsub("!","", column, fixed = TRUE))
  if(any(is.na(result) != is.na(column))) {
    stop("Problem converting")
  }
  result
}

functional_group_for_gene <- function(gene) {
  case_when(
    gene %in% sprintf("BBS%02d",c(1,2,4,5,7,8,9,18)) ~ "BBSome",
    gene == "BBS03" ~ "BBS03",
    gene %in% sprintf("BBS%02d",c(6,10,12)) ~ "Chaperonins",
    TRUE ~ "Others"
  )
}

age_transform_from_age <- function(age) {
  mean_age <- mean(age, na.rm=TRUE)
  sd_age <- sd(age, na.rm = TRUE)  
  function(x) {
    (x - mean_age)/ sd_age
  }
}

phenotype_long_map <- c(RD = "Retinal dystrophy",
                     OBE = "Obesity",
                     PD = "Polydactyly",
                     CI = "Cognitive impairment",
                     REP = "Reproductive system",
                     REN = "Renal anomalies",
                     HEART = "Heart anomalies",
                     LIV = "Liver anomalies",
                     DD = "Developmental delay")
  

phenotype_long_from_phenotype <- function(phenotype) {
  factor(phenotype_long_map[as.character(phenotype)], levels = phenotype_long_map[phenotypes_to_use])
}

read_data_base <- function(filename, sheet, ...) {
  read_excel(filename, sheet = sheet) %>%
    rename(case_no = "source case n.", 
           additional_mutations = "additional mutations", 
           mutation_types = "mut/mut",
           ethnic_group = "ethnic group",
           family_id = FamilyID) %>%
    filter(!is.na(source) | !is.na(gene)) %>% #NA in source is only in the empty rows at the end of the table
    select(source, case_no, family_id, gene, mutation_types, 
           sex, age, ethnic_group, ethnicity, ... , RD:DD) %>%
    mutate(
      CI = convert_uncertain_phenotype(CI),
      LIV = convert_uncertain_phenotype(LIV),
      REN = convert_uncertain_phenotype(REN),
      REP = convert_uncertain_phenotype(REP),
      DD = convert_uncertain_phenotype(DD),
      family_id = factor(family_id, exclude = "NA"),
      sex = factor(toupper(sex)), 
      source = factor(source),
      loss_of_function = factor(case_when(
        is.na(mutation_types) ~ "unknown",
        mutation_types == "trunc/trunc" ~ "certain",
        TRUE ~ "unknown"
      ), levels = c("unknown","certain")),
      loss_of_function_certain = as.numeric(loss_of_function == "certain")
    ) %>%
    rowid_to_column("ID") %>%
    #Code BBS to help ordering
    mutate(gene = factor(gsub("BBS([0-9])$","BBS0\\1", gene))) %>%
    
    
    #Introduce functional groups
    mutate(functional_group = 
             functional_group_for_gene(gene) %>% factor()
    ) 
  
}

read_main_data <- function() {
  
    data <- read_data_base(here("data","UPDATE8 SuppInfo Table S3 dataset-obn.xlsx"), sheet = "List1") %>%  
    
    #Get age categories
    mutate(age_corr = if_else(age == "5 month", as.character(5/12), gsub(",",".", age)),
           age_numbers = if_else(grepl("^[0-9]*\\.?[0-9]*$",age_corr), age_corr, NA_character_) %>% as.numeric(),
           age_group = factor(case_when(
             is.na(age_corr) ~ NA_character_,
             !is.na(age_numbers) ~ 
               case_when(
                 age_numbers < 10 ~ "0-9",
                 age_numbers < 20 ~ "10-19",
                 age_numbers < 30 ~ "20-29",
                 age_numbers < 40 ~ "30-39",
                 age_numbers < 50 ~ "40-49",
                 age_numbers < 60 ~ "50-59",
                 TRUE ~ "60+",
               ),
             TRUE ~ age_corr
           ))
    ) %>%
    
    #Get age for simpler model where categories are translated
    mutate(age_numbers_groups_guessed = 
             case_when(
               !is.na(age_numbers) ~ age_numbers,
               is.na(age_corr) ~ NA_real_,
               age_corr == "0-9" ~ 5,
               age_corr == "10-19" ~ 15,
               age_corr == "20-29" ~ 25,
               age_corr == "30-39" ~ 35,
               age_corr == "40-49" ~ 45,
               age_corr == "50-59" ~ 55,
               age_corr == "60+" ~ 70
             ),
           age_std_for_model = age_transform_from_age(age_numbers_groups_guessed)(age_numbers_groups_guessed)
    ) 
    
    
  
  if(!all(is.na(data$age) == is.na(data$age_group))) {
    stop("Error in age transforms")
  }
  if(!all(is.na(data$age) == is.na(data$age_numbers_groups_guessed))) {
    stop("Error in age transforms")
  }
  if(!all(is.na(data$age) == is.na(data$age_std_for_model))) {
    stop("Error in age transforms")
  }
  if(!(identical(suppressWarnings(as.numeric(data$age_corr)), data$age_numbers))) {
    stop("Error in age transforms")
  }
  
  data
}

read_validation_data <- function() {
  read_data_base(here("private_data","2019-04-16 London data for revision v04.xlsx"), sheet = "Sheet1", eGFR)
}

data_long_from_data  <- function(data) {
  data_long_all <- data %>% 
    gather("phenotype","phenotype_value",RD:DD)

  functional_groups_to_show <- data$functional_group %>% unique()
  
  data_long <- data_long_all %>%
    filter(phenotype %in% phenotypes_to_use) %>%
    filter(!is.na(phenotype_value)) %>%
    mutate(phenotype_value = as.integer(if_else(phenotype_value == 0, 0 , 1)),
           phenotype_long = phenotype_long_from_phenotype(phenotype),
           phenotype = factor(phenotype, levels = phenotypes_to_use)
           )
  
  if(any(is.na(data_long$phenotype))) {
    stop("Some phenotypes are NA")
  }
    
    
  data_long  
}