

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


read_main_data <- function() {
  
  data <- read_excel(here("private_data","EV table 2 dataset new.xlsx"), sheet = "List1") %>%
    rename(case_no = "source case n.", additional_mutations = "additional mutations", mutation_types = "mut/mut") %>%
    filter(!is.na(source) | !is.na(gene)) %>% #NA in source is only in the empty rows at the end of the table
    select(source, case_no,gene, mutation_types, 
           sex, age, RD:LIV) %>%
    mutate(
           CI = convert_uncertain_phenotype(CI),
           LIV = convert_uncertain_phenotype(LIV),
           REN = convert_uncertain_phenotype(REN),
           REP = convert_uncertain_phenotype(REP),
           Sex = factor(toupper(sex)), 
           source = factor(source),
           loss_of_function = factor(case_when(
                                     is.na(mutation_types) ~ "unknown",
                                     mutation_types == "trunc/trunc" ~ "certain",
                                     TRUE ~ "unknown"
                                     ), levels = c("unknown","certain")),
           loss_of_function_certain = as.numeric(loss_of_function == "certain")
    ) %>%
    rowid_to_column("ID") %>%  
    
    #Get age categories
    mutate(age_corr = if_else(age == "5 month", as.character(5/12), gsub(",",".", age)),
           age_numbers = if_else(grepl("^[0-9]*\\.?[0-9]*$",age_corr), age_corr, NA_character_) %>% as.numeric(),
           age_group = factor(case_when(
             is.na(age_corr) ~ NA_character_,
             !is.na(age_numbers) ~ 
               case_when(
                 age_numbers < 10 ~ "0-9",
                 age_numbers < 20 ~ "10-19",
                 age_numbers < 40 ~ "20-39",
                 TRUE ~ "40+",
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
               age_corr == "20-39" ~ 30,
               age_corr == "40+" ~ 50
             ),
           age_std_for_model = (age_numbers_groups_guessed - mean(age_numbers_groups_guessed, na.rm=TRUE))/ sd(age_numbers_groups_guessed, na.rm = TRUE)
    ) %>%
    
    
    #Code BBS to help ordering
    mutate(gene = factor(gsub("BBS([0-9])$","BBS0\\1", gene))) %>%
    
    
    #Introduce functional groups
    mutate(functional_group = 
             functional_group_for_gene(gene) %>% factor()
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


data_long_from_data  <- function(data) {
  data_long_all <- data %>% 
    gather("phenotype","phenotype_value",RD:LIV)

  functional_groups_to_show <- data$functional_group %>% unique()
  
  data_long <- data_long_all %>%
    filter(phenotype %in% phenotypes_to_use) %>%
    filter(!is.na(phenotype_value)) %>%
    mutate(phenotype_value = as.integer(if_else(phenotype_value == 0, 0 , 1)),
           phenotype = factor(phenotype, levels = phenotypes_to_use)
           )
  
  if(any(is.na(data_long$phenotype))) {
    stop("Some phenotypes are NA")
  }
    
    
  data_long  
}