.cran_packages = c("tidyverse", "psych", "cowplot", "DescTools")

## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

source("code/src/rfunctions.R")
source("code/src/ggstyles.R")
theme_set(mytheme(base_size = 12))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bridging.obj = read_rds("data/bridging_obj.RDS")

sc.pl = read.csv(
  "data/cpc_data.csv"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S1 : 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S : Numbers of circulating plasma cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sc.pl = sc.pl %>% 
  mutate(BT = factor(case_when(
    THERAPY_PRIOR_APHERESIS=="Talquetamab" ~ "Talquetamab",
    THERAPY_PRIOR_APHERESIS=="Teclistamab" ~ "Teclistamab",
    .default = "Other BT"), levels = c("Teclistamab", "Talquetamab", "Other BT")
  ))

lab = as_labeller(c("Pre Bite" = "pre BsAb", "LA" = "LA", "Day 30" = "Day 30", "Day 100" = "Day 100"))

numbers_of_cpc = sc.pl %>%
  mutate(Day = factor(Day, levels = c("Pre Bite", "LA", "Day 30", "Day 100"))) %>%
  group_by(Day, BT) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  ggplot(aes(x=BT, y=n, fill=BT)) +
  geom_bar(stat="identity", alpha=.8, color="black") +
  facet_wrap(~Day, nrow=1, labeller = lab) +
  mytheme_grid(12) +
  ylab("# cPC") +
  scale_fill_manual(name="Bridging", values=c("#994455", "#D93F00", "grey60")) +
  theme(axis.title.x = element_blank(), panel.grid.major.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1, vjust=1))

ggsave2(
  "figures/supplement/numbers_of_cpc_bars.png",
  numbers_of_cpc,
  height=70, width=110, units = "mm", bg="white", dpi=400, scale=1.5
)
