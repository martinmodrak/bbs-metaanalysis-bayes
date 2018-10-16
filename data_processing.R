

convert_czech_numeric <- function(column) {
  result <- as.numeric(gsub(",",".", column, fixed = TRUE))
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
  
  data <- read_excel(here("private_data","--META-DATA-5.xlsx"), sheet = "DATA 5") %>%
    rename(case_no = "source case n.", mut_nucl = "mu nucl", mut_amino = "principal mutation", mut_type = "mu-type", additional_mutation = "additional mutation", REPROD = "Reprod. organs abnorm") %>%
    select(source, case_no,gene, mut_nucl, mut_amino, mut_type, 
           Sex, Age, RD:speech) %>%
    mutate(CI = convert_czech_numeric(CI), DD = convert_czech_numeric(DD),
           LIV = convert_czech_numeric(LIV), REN = convert_czech_numeric(REN),
           REPROD = convert_czech_numeric(REPROD),
           Sex = factor(toupper(Sex)), 
           source = factor(source)
    ) %>%
    rowid_to_column("ID") %>%  
    
    #Get age categories
    mutate(Age_corr = if_else(Age == "5 month", as.character(5/12), gsub(",",".", Age)),
           age_numbers = if_else(grepl("^[0-9]*\\.?[0-9]*$",Age_corr), Age_corr, NA_character_) %>% as.numeric(),
           age_group = factor(case_when(
             is.na(Age_corr) ~ NA_character_,
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
             TRUE ~ Age_corr
           ))
    ) %>%
    
    #Get age for simpler model where categories are translated
    mutate(age_numbers_groups_guessed = 
             case_when(
               !is.na(age_numbers) ~ age_numbers,
               is.na(Age_corr) ~ NA_real_,
               Age_corr == "0-9" ~ 5,
               Age_corr == "10-19" ~ 15,
               Age_corr == "20-29" ~ 25,
               Age_corr == "30-39" ~ 35,
               Age_corr == "40-49" ~ 45,
               Age_corr == "50-59" ~ 55,
               Age_corr == "60+" ~ 65
             ),
           age_std_for_model = (age_numbers_groups_guessed - mean(age_numbers_groups_guessed, na.rm=TRUE))/ sd(age_numbers_groups_guessed, na.rm = TRUE)
    ) %>%
    
    
    #Code BBS to help ordering
    mutate(gene = factor(gsub("BBS([0-9])$","BBS0\\1", gene))) %>%
    
    
    #Introduce functional groups
    mutate(functional_group = 
             functional_group_for_gene(gene) %>% factor()
    ) 
  
  if(!all(is.na(data$Age) == is.na(data$age_group))) {
    stop("Error in age transforms")
  }
  if(!all(is.na(data$Age) == is.na(data$age_numbers_groups_guessed))) {
    stop("Error in age transforms")
  }
  if(!all(is.na(data$Age) == is.na(data$age_std_for_model))) {
    stop("Error in age transforms")
  }
  if(!(identical(suppressWarnings(as.numeric(data$Age_corr)), data$age_numbers))) {
    stop("Error in age transforms")
  }
  
  #For some reasone this cannot be done with dplyr
  data$hyposmia_anosmia = data$`hyposmia/anosmia`
  data$`hyposmia/anosmia` = NULL
  
  data
}


data_long_from_data  <- function(data) {
  data_long_all <- data %>% 
    gather("phenotype","phenotype_value",RD:speech) 
  
  data_long_all %>% group_by(phenotype) %>% summarise(count = sum(!is.na(phenotype_value)), prop = mean(phenotype_value, na.rm = TRUE))
  
  
  functional_groups_to_show <- data$functional_group %>% unique()
  
  data_long <- data_long_all %>%
    filter(phenotype %in% phenotypes_to_use) %>%
    filter(!is.na(phenotype_value)) %>%
    mutate(phenotype_value = as.integer(if_else(phenotype_value == 0, 0 , 1)))
  
  data_long  
}