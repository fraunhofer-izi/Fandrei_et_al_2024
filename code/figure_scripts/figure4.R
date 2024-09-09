.cran_packages = c("tidyverse", "psych", "cowplot",
                    "ggpubr", "ggh4x", "ggpmisc", "DescTools")

## Loading library
for (pack in .cran_packages) {
  suppressMessages(library(
    pack,
    quietly = TRUE,
    verbose = FALSE,
    character.only = TRUE
  ))
}

source("code/src/ggstyles.R")
theme_set(mytheme(base_size = 12))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD DATA
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bridging.obj = read_rds("../data/bridging_obj.RDS")

clin_data = bridging.obj$clin_data
immune_df = bridging.obj$immune_df
immune_subset_df = bridging.obj$immune_subsets
cytotox_df = bridging.obj$cytotox
elisa_df = bridging.obj$elisa_df


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 4A : CAR expansion stratified by bridging
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

car_expansion_bridging <- immune_subset_df %>%
  filter(!Day %in% c("Leukapheresis", "Day 0"), immune_subset %in% c("CD3_CAR_%", "CD4_CAR_%", "CD8_CAR_%")) %>%
  mutate(immune_subset = car::recode(immune_subset, recode_car_subsets)) %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill=therapy_before_cart)) +
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
  scale_fill_manual("Bridging", values=anno_nejm) +
  scale_y_continuous(breaks=c(0,25,50,75,100), expand=c(0.15,0)) +
  stat_compare_means(label="p", label.y=120, hide.ns=T, size=3.5, hjust=0.3, vjust=-.4) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, label.size=5.5)

ggsave2(
  "figures/main/figure_4/car_expansion_bridging.png",
  car_expansion_bridging,
  width = 120, height = 100, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# CD4 CAR / CD8 CAR ratio boxplots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Calculate ratios for CAR T and total Total T cells
recode_car_subsets <- "'CD3_CAR_%' = 'CD3+CAR+\\n(% of CD3+)'; 'CD4_CAR_%'='CD4+CAR+\\n(% of CD3+CAR+)'; 'CD8_CAR_%' ='CD8+CAR+\\n(% of CD3+CAR+)'"

ratio_df <- immune_df %>% 
  mutate(across(c(`CD4_CAR_%`, `CD8_CAR_%`), function(x) gsub("n.d.", "", x))) %>%
  mutate(across(c(`CD4_CAR_%`, `CD8_CAR_%`), as.numeric)) %>%
  mutate(ratio = `CD4_CAR_%` / `CD8_CAR_%`) %>%
  dplyr::select(patient_id, Day, ratio) %>%
  filter(!Day %in% c("Leukapheresis", "Day 0")) %>%
  left_join(clin_data[,c("patient_id", "30 day response", "status_before_cart", "therapy_before_cart", "Product", "refractoriness_group")], by="patient_id") %>%
  mutate(bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other"), bridging_num = ifelse(therapy_before_cart == "Bispecific Ab", 1, 0))

recode_total_cd3 <- "'CD4_abs' = 'CD4'; 'CD8_abs' = 'CD8'"

ratio_total_cd3 <- immune_df %>% 
  mutate(across(c(`CD4_abs`, `CD8_abs`), function(x) gsub("n.d.", "", x))) %>%
  mutate(across(c(`CD4_abs`, `CD8_abs`), as.numeric)) %>%
  mutate(ratio = `CD4_abs` / `CD8_abs`) %>%
  dplyr::select(patient_id, Day, ratio) %>%
  left_join(clin_data[,c("patient_id", "30 day response", "status_before_cart", "therapy_before_cart", "Product", "refractoriness_group")], by="patient_id") %>%
  mutate(bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other"))

# Boxplot
ratio_bridging_binary = ratio_df %>%
  ggplot(aes(x=Day, y=ratio, fill=bridging_bin)) +
  geom_boxplot(alpha=.75, outlier.shape = NA, show.legend = F) +
  geom_point(position=position_jitterdodge(jitter.width = 0.15), alpha=.7, show.legend = F, size=.9) +
  scale_fill_manual("Bridging", values=c(anno_nejm[1], Other="grey60")) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  ) +
  scale_y_continuous(expand=c(0,0.5)) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, method.args = list(alternative = "greater", paired = FALSE)) +
  ylab("CD4+CAR+/CD8+CAR+ ratio") +
  guides(fill=guide_legend(nrow = 1))

tapply(ratio_df[ratio_df$Day=="Day 7",]$ratio, ratio_df[ratio_df$Day == "Day 7",]$bridging_bin, summary)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 4B : Cytotox per bridging
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cytotox_plot <- cytotox_df %>%
  mutate(condition = car::recode(condition, "'U266+T'='U266 + T'; 'U266+CAR-T'='U266 + CAR-T'; 'U266'='U266 control'")) %>%
  mutate(condition=factor(condition, levels=c("U266 control", "U266 + T", "U266 + CAR-T"))) %>%
  filter(timepoint=="24h") %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill=therapy_before_cart)) +
  geom_boxplot(alpha=.75, outlier.shape = NA) +
  geom_point(color="grey10", position = position_jitterdodge(jitter.width=0.4), alpha=.7, size=1.3) + 
  facet_grid(.~condition) +
  mytheme_grid(13) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual("Bridging", values=anno_nejm) +
  scale_y_continuous(breaks=seq(0,100,25), expand=c(0.1,0)) +
  ylab("% Viable Cells") +
  guides(fill=guide_legend(nrow=1)) +
  stat_compare_means(label = "p.format", vjust=-4) +
  geom_pwc(label="p.signif", hide.ns="p", vjust=.5)

cytotox_df %>%
  filter(timepoint=="24h") %>%
  group_by(therapy_before_cart, condition) %>%
  summarize(median=median(perc, na.rm=T), min=min(perc, na.rm=T), max=max(perc, na.rm=T))

ggsave2(
  "figures/main/figure_4/cytotox_bridging.png",
  height=80, width=98, units="mm", dpi=400, bg="white", scale=1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 4C : CD4 / CD8  ratio longitudinal
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ratio_df %>%
  group_by(Day, bridging_bin) %>%
  summarize(median = median(ratio, na.rm=T), min = min(ratio, na.rm=T), max = max(ratio, na.rm=T))

# CAR
longitudinal_cart_ratio_bite_p <- JonckheereTerpstraTest(Day ~ ratio, subset(ratio_df, subset = bridging_bin == "Bispecific Ab"))
longitudinal_cart_ratio_other_p <- JonckheereTerpstraTest(Day ~ ratio, subset(ratio_df, subset = bridging_bin == "Other"))

cart_jockheere_df <- data.frame(
  bridging_bin = c("Bispecific Ab", "Other"),
  p_signif =  c(round(longitudinal_cart_ratio_bite_p$p.value, 3), round(longitudinal_cart_ratio_other_p$p.value, 3))
)

violins_cart <- ratio_df %>%
  ggplot(aes(x=Day, y=ratio, fill=bridging_bin)) +
  # geom_boxplot(alpha=.75, outlier.shape = NA) +
  geom_violin(draw_quantiles = .5, alpha = .75, color = "grey20", linetype="longdash") +
  geom_point(color="grey10", position = position_jitter(width=0.3), alpha=.7, size = .8, show.legend = F) +
  scale_fill_manual("Bridging", values=c(anno_nejm[1], Other = "grey60")) +
  facet_wrap(.~bridging_bin, nrow=1) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle=45, hjust=1)
  ) +
  geom_bracket(
    data = cart_jockheere_df, aes(label = paste("P=", p_signif)),
    xmin = "Day 7", xmax = "Day 100", y.position = 4,
  ) +
  scale_y_continuous(expand=c(0,0.5)) +
  stat_compare_means(label="p", hide.ns=T, size=4.5) +
  xlab("Time point") +
  ylab("CD4+CAR+/CD8+CAR+ ratio") +
  guides(fill=guide_legend(nrow = 1)) +
  geom_hline(yintercept=1, alpha=.4, linetype = "dotted")

options(scipen=0)
# Total CD3
longitudinal_ratio_bite_p <- JonckheereTerpstraTest(Day ~ ratio, subset(ratio_total_cd3, subset = bridging_bin == "Bispecific Ab"))
longitudinal_ratio_other_p <- JonckheereTerpstraTest(Day ~ ratio, subset(ratio_total_cd3, subset = bridging_bin == "Other"))
total_cd3_jockheere_df <- data.frame(
  bridging_bin = c("Bispecific Ab", "Other"),
  p_signif =  c(paste0("=", round(longitudinal_ratio_bite_p$p.value, 2)), ifelse(longitudinal_ratio_other_p$p.value<0.001, "<0.001", signif(longitudinal_ratio_other_p$p.value, digits=2)))
)

violins_total_cd3 <- ratio_total_cd3 %>%
  mutate(Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                      levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  ggplot(aes(x=Day, y=ratio, fill=bridging_bin)) +
  geom_boxplot(alpha=.75, outlier.shape = NA) +
  geom_point(color="grey10", position = position_jitter(width=0.3), alpha=.7, size = .8, show.legend = F) +
  scale_fill_manual("Bridging", values=c(anno_nejm[1], Other = "grey60")) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45, hjust=1)
  ) +
  scale_y_continuous(expand=c(0.1,0), trans="log10") +
  stat_compare_means(label="p", hide.ns="p", size=4.5) +
  xlab("Time point") +
  ylab("Total CD4/CD8 ratio") +
  geom_pwc(label="p.format", hide.ns="p") +
  guides(fill=guide_legend(nrow = 1))
violins_total_cd3

ratio_violins <- ratio_bridging_binary + violins_total_cd3 +
  plot_layout(guides = "collect",  nrow=2) & theme(legend.position = "bottom")
ratio_violins

ggsave2(
  filename="figures/main/figure_4/ratio_longitudinal.png",
  ratio_violins,
  width = 115, height = 160, dpi = 400, bg = "white", units = "mm", scale = 1.1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 4D: serum BCMA levels
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

elisa_df[,c("patient_id", "Day", "BCMA [ng/mL]")] %>%
  mutate(sbcma = `BCMA [ng/mL]`) %>%
  left_join(clin_data[,c("patient_id", "treatment_before_apheresis", "RESPONSE_CONSENSUS_1")]) %>%
  filter(!Day %in% c("Day 7" , "Day 14")) %>%
  mutate(group = factor(case_when(
    treatment_before_apheresis == "Teclistamab" ~ "teclistamab",
    treatment_before_apheresis == "Talquetamab" ~ "talquetamab",
    .default = "other"
  ), levels=c("talquetamab","teclistamab", "other")) ) %>%
  group_by(group, Day) %>%
  summarize(median=median(sbcma, na.rm=T), min=min(sbcma, na.rm=T), max=max(sbcma, na.rm=T)) %>%
  filter(Day == "Day 0")

elisa_df[,c("patient_id", "Day", "BCMA [ng/mL]")] %>%
  left_join(clin_data[,c("patient_id", "treatment_before_apheresis", "RESPONSE_CONSENSUS_1")]) %>%
  filter(!Day %in% c("Day 7" , "Day 14")) %>%
  mutate(group = factor(case_when(
    treatment_before_apheresis == "Teclistamab" ~ "teclistamab",
    treatment_before_apheresis == "Talquetamab" ~ "talquetamab",
    .default = "other"
  ), levels=c("talquetamab","teclistamab", "other")) ) %>%
  ggplot(aes(x=Day, y=`BCMA [ng/mL]`+0.1, fill=group)) +
  geom_boxplot(outlier.shape = NA, alpha=.75) +
  geom_point(position = position_jitterdodge(jitter.width = .15), alpha=.5, size=.9, outlier.shape=NA) +
  scale_fill_manual(name="Bridging", values=c("#EE7733", "#994455", "grey60")) +
  scale_y_continuous(expand=c(0.15,0), trans="log10") +
  stat_compare_means(label="p.format", vjust=-2.7) +
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
  "figures/main/figure_4/sbcma_bridging.png",
  width = 95, height = 84, dpi = 400, bg = "white", units = "mm", scale = 1.6
)
