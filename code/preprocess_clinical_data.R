library(tidyverse)
library(readxl)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Import Data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Main clinical dataframe
clin_data <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=1,
  na="NA"
)

# PCR data
pcr <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=3,
  skip=1
)

# Refractoriness
refractoriness_df <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=5
)

# Cytotoxicity
cytotox_data <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=8,
  na = "NA"
)

# Immune subsets
immune_status <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=6
)

# T cell subsets
t_cells <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=7
)

# CRS
crs <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=4,
  na = "NA"
)

# BITE administration
bite_admin_df <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=9,
  na = "NA"
)

# BITE leukapheresis
bite_leukapheresis <- read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=10
)

# Cytotoxicity
cytotoxicity_immu_df = read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=11
)

# Immune checkpoints FC
checkpoint_immu_df = read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=12
)

# PD1 CARs
pd1_car_df = read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=13
)

# ELISA
elisa_df = read_excel(
  "data/cart_bridging_data.xlsx",
  sheet=14
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Factors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Day levels
day_levels = c("Leukapheresis", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
# Response groups
response_groups_str <- "'nCR'='CR'; 'VGPR'='VGPR/PR'; 'PR'='VGPR/PR'; 'SD'='SD/PD'; 'PD'='SD/PD'"
# Response group levels
response_levels = c("CR", "VGPR/PR", "SD/PD")
# Bridging therapy levels
therapy_levels = c('Bispecific Ab', 'Chemotherapy', 'CD38', 'SLAMF7')

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Preprocess
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Clinical data
clin_data <- clin_data %>%
  mutate(across(c(13:33), ~ifelse(. == "x", 1, 0)))  %>%
  mutate(
    patient = gsub("atient0", "", patient_id),
    therapy_before_cart = factor(therapy_before_cart),
    Product = car::recode(type_cart, "'a'='Ide-cel'; 'c'='Cilta-cel'", as.factor=T, levels=c("Ide-cel", "Cilta-cel")),
    sex = factor(sex),
    crs_grade = car::recode(crs_grade, "'0'='no CRS'; '1'='Grade 1'; '2'='Grade 2'", as.factor = T, levels = c("no CRS", "Grade 1", "Grade 2")),
    time_diag_cart = time_diag_cart/30.417,
    # convert cytogenetics to binary
    mutate(across(c(10:11, 13:33), function(x) {ifelse(is.na(x), 0, x)}))
  )
  
clin_data <- clin_data %>% 
  mutate(
    sample_id = paste0(patient_id, "_1"),
    `30 day binary` = factor(ifelse(`30 day response` %in% c("nCR", "CR", "VGPR", "PR"), "Responder", "Non-Responder"), levels = c("Responder", "Non-Responder")),
    `30 day response` = car::recode(`30 day response`, response_groups_str),
    `3 month response` = car::recode(`3 month response`, response_groups_str),
    all_groups_before_cart = status_before_cart,
    status_before_cart = car::recode(status_before_cart, response_groups_str),
    riss = as.factor(`R-ISS`),
    # convert time columns to months
    os_months = os/30.417,
    pfs_months=pfs/30.417,
    bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other")
  )

clin_data$`30 day response` <- factor(clin_data$`30 day response`, levels = response_levels)
clin_data$`3 month response` <- factor(clin_data$`3 month response`, levels = response_levels)
clin_data$RESPONSE_CONSENSUS_1 = ifelse(clin_data$`30 day response`=="CR", "CR", "non-CR")
clin_data$status_before_cart <- factor(clin_data$status_before_cart, levels = response_levels)
clin_data$therapy_before_cart <- factor(clin_data$therapy_before_cart, levels=therapy_levels)

# Add prior lines of treatment
lines_treatment_df <- clin_data %>% 
  select(c(1, grep("Regime", colnames(clin_data)))) %>%
  pivot_longer(!patient_id, names_to="regime", values_to="agent") %>%
  filter(!is.na(agent)) %>%
  group_by(patient_id) %>%
  filter(row_number() != n()) %>%
  summarize(n_lines_treatment=n())

clin_data = left_join(clin_data, lines_treatment_df, by="patient_id")

# PCR data
pcr_data <- pcr[,c(1, 5:9)] %>%
  pivot_longer(!patient_id, names_to = "sample", values_to = "no_copies")

pcr_data$Day <- lapply(pcr_data$sample, function(x) last(str_split_1(x, "_"))) %>% unlist()

pcr_data <- pcr_data %>%
  mutate(Day = factor(paste0("Day ", Day), levels=c("Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  left_join(clin_data[,c("patient_id", "status_before_cart", "therapy_before_cart", "Product", "30 day response",
                         "refractoriness_group")], by="patient_id")

# Refractoriness
refractoriness_df <- refractoriness_df %>%
  mutate(refractoriness_group = case_when(
    penta == "yes" ~ "PentaRRMM",
    `triple class` == "refractory" ~ "TCRRMM",
    `triple class` == "exposed" ~  "TCexposed"),
    patient_id = Patient
  )

# Cytotoxicity
cytotox_df <- cytotox_data %>%
  pivot_longer(c(2:5)) %>%
  left_join(clin_data[,c("patient_id", "therapy_before_cart", "status_before_cart", "Product", "30 day response")], by="patient_id")

# Immune subsets
immune_df <- immune_status %>%
  left_join(clin_data[,c("patient_id", "30 day response",
                         "status_before_cart", "therapy_before_cart", "Product", "refractoriness_group")], by="patient_id") %>%
  filter(!is.na(patient_id)) %>%
  group_by(patient_id) %>%
  mutate(
    sample_id = paste0(patient_id, "_", row_number()),
    Day = factor(Day, levels=day_levels),
    `30 day response` = factor(car::recode(`30 day response`, response_groups_str), levels = response_levels),
    status_before_cart = factor(car::recode(status_before_cart, response_groups_str), levels = response_levels)
  ) %>%
  ungroup()

## Columns with % values for subsets
immune_perc_col <- grep("?%", colnames(immune_status), value=T)

immune_subset_df <- immune_df %>%
  dplyr::select(all_of(c("sample_id", immune_perc_col))) %>%
  mutate(across(c(2:7), as.numeric)) %>%
  pivot_longer(!sample_id, names_to = "immune_subset", values_to = "perc") %>%
  rowwise() %>%
  mutate(patient_id = str_split_1(sample_id , "_")[1]) %>%
  ungroup() %>%
  left_join(immune_df[,c("sample_id", "Day", "30 day response", "status_before_cart", "therapy_before_cart",
                         "Product", "refractoriness_group")], by="sample_id")

# T cell subsets
recode_tcells <- "'CD4 Thymusemi %' = 'CD4+ thymic emigrants'; 'naive CD4 %' = 'Naive CD4+'; 'CD4 Eff-Memory %' = 'CD4+ effector memory';
    'CD4 central Memory %' = 'CD4+ central memory'; 'CD4 Effektorzellen %' = 'CD4+ effector'; 'CD8 Thymusemigrants %' = 'CD8+ thymic emigrants';
    'naive CD8 %' = 'Naive CD8+'; 'CD8 Eff-Memory %' ='CD8+ effector memory'; 'CD8 central Memory %' = 'CD8+ central memory';
    'CD8 Effektorzellen %' = 'CD8+ effector'; 'Tregs %' = 'Tregs';  'naive Tregs %' = 'Naive Tregs';
    'Memory Tregs %'= 'Memory Tregs'"

t_cells$Day = factor(t_cells$Day, levels = day_levels)

t_cells_long = t_cells %>%
  pivot_longer(cols = c(4:ncol(t_cells)), names_to = "tcell_subset") %>%
  mutate(group = case_when(
    grepl("CD4", tcell_subset) ~ "CD4",
    grepl("CD8", tcell_subset) ~ "CD8",
    grepl("Treg", tcell_subset) ~ "Treg"
  ))

t_cells_abs = t_cells_long %>%
  filter(grepl("abs", tcell_subset)) %>%
  mutate(tcell_subset = car::recode(gsub("abs", "%", tcell_subset), recode_tcells) ) %>%
  left_join(clin_data[,c("patient_id", "therapy_before_cart", "status_before_cart", "bridging_bin")])

# CRS
crs_df <- crs %>% 
  mutate(crs_group = factor(crs_group, levels=c("no", "yes-Toci", "yes+Toci")), patient_id = gsub("_1", "", sample_id)) %>%
  left_join(clin_data[,c("sample_id", "status_before_cart", "therapy_before_cart")], by="sample_id")

# Checkpoints
checkpoint_immu_df = checkpoint_immu_df %>%
  filter(!Day %in% c("Lymphodepletion", "3")) %>%
  mutate(Day = factor(
    ifelse(Day == "Leukapheresis", "LA", paste0("Day ", Day)), levels = c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
  ))
 
# # Cytotoxicity Immu
cytotoxicity_immu_df = cytotoxicity_immu_df %>%
  mutate(across(3:17, as.numeric)) %>%
  filter(!Day %in% c("Lymphodepletion", "3")) %>%
  mutate(Day = factor(
    ifelse(Day=="Leukapheresis", "LA", paste0("Day ", Day)),
           levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
  )) %>%
  left_join(clin_data[,c("patient_id", "bridging_bin")], by = "patient_id") %>%
  select(patient_id, bridging_bin, Day, (grep("CD8+", colnames(cytotoxicity_immu_df))), (grep("CD4+", colnames(cytotoxicity_immu_df)))) %>%
  pivot_longer(cols = 4:15) %>%
  mutate(
    CAR = ifelse(gsub(".+CAR", "", name) == "+", "CAR+", "CAR-"),
    marker = gsub("\\+CAR[\\+,\\-]", "", gsub("CD[4,8]\\+", "", name)),
    line = str_sub(name, start=1, end=4)
  )

# Elisa
elisa_df = elisa_df %>%
  mutate(Day = factor(
    case_when(
      Day=="Leukapheresis" ~ "LA", 
      Day == "pre_BT" ~ "pre BT",
      .default=Day
    ),
  levels = c("pre BT", "LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Compose RDS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bridging.obj <- list(
  clin_data = clin_data,
  pcr_data = pcr_data,
  refractoriness = refractoriness_df,
  immune_df = immune_df,
  immune_subsets = immune_subset_df,
  cytotox = cytotox_df,
  t_cells = t_cells,
  t_cells_abs = t_cells_abs,
  crs = crs_df,
  bite_admin = bite_admin_df,
  bite_leukapheresis = bite_leukapheresis,
  checkpoint_immu_df = checkpoint_immu_df,
  pd1_car_df = pd1_car_df,
  cytotoxicity_immu_df = cytotoxicity_immu_df,
  elisa_df = elisa_df
)

write_rds(bridging.obj, "data/bridging_obj.RDS")
