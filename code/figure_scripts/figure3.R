.cran_packages = c("tidyverse", "psych", "cowplot", "ggpubr", "ggh4x", "ggpmisc")

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

bridging.obj = read_rds("../data/bridging_obj.RDS")

clin_data <- bridging.obj$clin_data
immune_df = bridging.obj$immune_df
immune_subset_df = bridging.obj$immune_subsets
cytotoxicity_immu_df = bridging.obj$cytotoxicity_immu_df
elisa_df = bridging.obj$elisa_df

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CD4 CAR / CD8 CAR ratio
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

recode_car_subsets <- "'CD3_CAR_%' = 'CD3+CAR+\\n(% of CD3+)'; 'CD4_CAR_%'='CD4+CAR+\\n(% of CD3+CAR+)'; 'CD8_CAR_%' ='CD8+CAR+\\n(% of CD3+CAR+)'"

ratio_df <- immune_df %>% 
  mutate(across(c(`CD4_CAR_%`, `CD8_CAR_%`), function(x) gsub("n.d.", "", x))) %>%
  mutate(across(c(`CD4_CAR_%`, `CD8_CAR_%`), as.numeric)) %>%
  mutate(ratio = `CD4_CAR_%` / `CD8_CAR_%`) %>%
  dplyr::select(patient_id, Day, ratio) %>%
  filter(!Day %in% c("Leukapheresis", "Day 0")) %>%
  left_join(clin_data[,c("patient_id", "30 day response", "status_before_cart", "therapy_before_cart", "Product", "refractoriness_group")], by="patient_id")

ratio_df %>%
  mutate(response=`30 day response`) %>%
  group_by(Day, response) %>%
  summarize(median = median(ratio, na.rm=T), min = min(ratio, na.rm=T), max = max(ratio, na.rm=T))

ratio_30_day_response = ratio_df %>%
  ggplot(aes(x=`30 day response`, y=ratio, fill=`30 day response`)) +
  geom_boxplot(alpha=.75, outlier.shape = NA) +
  geom_point(color="grey10", position = position_jitter(width=0.3), alpha=.7, size=1) +
  scale_fill_manual(values=palette_response) +
  facet_wrap(.~Day, nrow=1) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    # panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
  ) +
  scale_y_continuous(expand=c(0,0.5)) +
  stat_compare_means(label="p", hide.ns=T, size=4.5, label.y = 4.7) +
  geom_pwc(label="p.signif", hide.ns="p") +
  ylab("CD4+CAR+/CD8+CAR+ ratio") +
  guides(fill=guide_legend(ncol=3, title.position = "left"))

ggsave2(
  "figures/main/figure_3/ratio_30_day_response.png",
  ratio_30_day_response,
  width = 120, height = 75, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3B : CAR expansion 30 day response
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

car_expansion_30_day_response = immune_subset_df %>%
  filter(!Day %in% c("Leukapheresis", "Day 0"), immune_subset %in% c("CD3_CAR_%", "CD4_CAR_%", "CD8_CAR_%")) %>%
  mutate(immune_subset = car::recode(immune_subset, recode_car_subsets)) %>%
  ggplot(aes(x=`30 day response`, y=perc, fill=`30 day response`)) +
  geom_boxplot(alpha=.8, outlier.shape = NA) +
  mytheme_grid(13) +
  geom_point(color="grey10", position = position_jitter(width=0.3), alpha=.7, size=1.2) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank()
  ) +
  facet_grid(immune_subset~Day, scales="free_y") +
  ylab("Proportion (%)") +
  scale_fill_manual(values=palette_response) +
  scale_y_continuous(breaks=c(0,25,50,75,100), expand=c(0.15,0)) +
  stat_compare_means(label="p", label.y=120, hide.ns=T, size=3.5, hjust=0.3, vjust=-.4) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, label.size=5.5)

ggsave2(
  "figures/main/figure_3/car_expansion_30_day_response.png",
  car_expansion_30_day_response,
  width = 120, height = 100, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3C : cytotoxicity markers on T cells
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cytotoxicity_immu_df %>%
  filter(!is.na(value), CAR=="CAR+", !Day %in% c("LA", "Day 0"), !marker %in% c("CAR+", "CAR-")) %>%
  left_join(clin_data[,c("patient_id", "RESPONSE_CONSENSUS_1")]) %>%
  mutate(marker = ifelse(marker=="Grz", "Granzyme", "Perforin")) %>%
  ggplot(aes(x=Day, y=value, fill=RESPONSE_CONSENSUS_1)) +
  geom_boxplot(outlier.shape=NA, alpha=.6) +
  geom_point(color="grey10", position = position_jitterdodge(jitter.width=0.3), alpha=.9, size=1.2) +
  facet_grid(line~marker) +
  mytheme_grid(13) +
  theme(
    axis.text.x= element_text(angle=45, hjust=1, vjust=1),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  ) +
  xlab("Time point after CAR T infusion") +
  ylab("Proportion (%) of CAR+ cells") +
  scale_fill_manual(name="", values=c("#6699CC", "#997700")) +
  scale_y_continuous(expand=c(0.1,0)) +
  geom_pwc(label="p.signif", hide.ns="p")

ggsave2(
  "figures/main/figure_3/cytotox_markers.png",
  width = 73, height = 90, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 3D : serum BCMA (sBCMA) levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

elisa_df[,c("patient_id", "Day", "BCMA [ng/mL]")] %>%
  left_join(clin_data[,c("patient_id", "bridging_bin", "RESPONSE_CONSENSUS_1")]) %>%
  mutate(bcma = `BCMA [ng/mL]`) %>%
  filter(Day != "Day 14") %>% 
  group_by(Day, RESPONSE_CONSENSUS_1) %>%
  summarize(median = median(bcma, na.rm=T), min = min(bcma, na.rm=T), max = max(bcma, na.rm=T))

elisa_df[,c("patient_id", "Day", "BCMA [ng/mL]")] %>%
  left_join(clin_data[,c("patient_id", "bridging_bin", "RESPONSE_CONSENSUS_1")]) %>%
  filter(Day != "Day 14") %>%
  ggplot(aes(x=Day, y=`BCMA [ng/mL]`+0.1, fill=RESPONSE_CONSENSUS_1)) +
  geom_boxplot(outlier.shape = NA, alpha=.75) +
  geom_point(position = position_jitterdodge(jitter.width = .15), alpha=.5, size=.9, outlier.shape=NA) +
  scale_fill_manual(name="", values=c("#6699CC", "#997700")) +
  scale_y_continuous(expand=c(0.1,0), trans="log10") +
  geom_pwc(label="p.signif", hide.ns="p") +
  ylab("sBCMA (ng/mL)") +
  xlab("Time point") +
  facet_grid(.~paste0("sBCMA")) +
  mytheme_grid(13) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave2(
  "figures/main/figure_3/sbcma_response_consensus.png",
  width = 80, height = 75, dpi = 400, bg = "white", units = "mm", scale = 1.5
)
