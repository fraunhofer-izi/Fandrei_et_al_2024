.cran_packages = c("tidyverse", "table1")


## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bridging.obj = read_rds("../data/bridging_obj.RDS")
clin_data =  bridging.obj$clin_data

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Statistical tests
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pvalue <- function(x, ...) {
  x <- x[-length(x)]  # Remove "overall" group
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform an ANOVA
    # p <- summary(aov(y ~ g))[[1]][["Pr(>F)"]][1]
    # For numeric variables, perform Wilcoxon's rank-sum test
    p <- wilcox.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           # "Median (SD)"=sprintf("%s (&plusmn; %s)", MEDIAN, SD),
                                                           "Median (IQR)"=sprintf("%s [%s - %s]", MEDIAN, Q1, Q3),
                                                           "Range" = sprintf("%s - %s", MIN, MAX)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Table 1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_data$gain_1q = as.factor(ifelse(clin_data$gain_1q | clin_data$gain_1q == 1, "1", "0"))
clin_data$del_17p = as.factor(clin_data$del_17p)
clin_data$`t(4;14)` = as.character(clin_data$`t(4;14)`)
clin_data$death = as.character(clin_data$death)
clin_data$refractoriness_group = factor(clin_data$refractoriness_group, levels = c("TCexposed", "TCRRMM", "PentaRRMM"))

labels <- list(
  variables=list(
    age_cart = "Age at CAR T therapy",
    sex="Sex",
    riss = "R-ISS risk category",
    `Salmon and Durie` = "Prognostic assessment according to Salmon and Durie",
    del_17p = "del_17p",
    `t(4;14)` = "t(4;14)",
    gain_1q = "gain1q",
    time_diag_cart = "Time from primary diagnosis to CAR T cell infusion (months)",
    time_apheresis_cart = "Time from apharesis to infusion (months)",
    status_before_cart = "Remission status prior to CAR T",
    therapy_before_cart = "Type of bridging therapy prior to CAR T",
    n_lines_treatment="Number of lines of treatment prior to CAR T",
    refractoriness_group = "Refractoriness prior to CAR T",
    `30 day response` = "Response on Day 30 after CAR T",
    crs_grade="CRS status post CAR T",
    death = "Death status"
  ),
  groups=list("", "Product")
)

strata_clin_data <- c(list(Total=clin_data), split(clin_data, clin_data$Product))

t1 <- table1(
  strata_clin_data,
  # groupspan=c(1, 2),
  labels = labels,
  render.continuous=my.render.cont,
  render.categorical=my.render.cat,
  extra.col = list('P-Value'=pvalue),
)
t1

write.csv(
  as.data.frame(t1),
  "figures/main/tables/table1.csv",
  quote=F,
  row.names=F
)
