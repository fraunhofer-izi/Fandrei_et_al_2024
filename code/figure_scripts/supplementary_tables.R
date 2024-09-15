library(tidyverse)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bridging.obj = read_rds("data/bridging_obj.RDS")

clin_data = bridging.obj$clin_data

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Table S1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

table_s1 = clin_data[,c("patient_id", "therapy_before_cart", "refractoriness_group", "status_before_cart", "30 day response", "pfs", "progression", "os", "death")]

write.csv(
  table_s1,
  "figures/supplement/table_s1.csv",
  quote=F,
  row.names=F
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Table S2
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

table_s2 = clin_data[,c("patient_id", "therapy_before_cart", "treatment_before_apheresis")]

write.csv(
  table_s2,
  "figures/supplement/table_s2.csv",
  quote=F,
  row.names=F
)
