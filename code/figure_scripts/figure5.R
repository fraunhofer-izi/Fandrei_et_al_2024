.cran_packages = c("tidyverse", "psych", "cowplot", "DescTools", "readxl")

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

clin_data = bridging.obj$clin_data
tcell.df = bridging.obj$t_cells
t_cells_abs_df = bridging.obj$t_cells_abs
pd1_car_df = bridging.obj$pd1_car_df
checkpoint_immu_df = bridging.obj$checkpoint_immu_df

# Single cell data
sc.t = read.csv(
  "data/tcell_data.csv"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#  Fig. 5A : Tregs in BT groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tregs_bridging = t_cells_abs_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))
  ) %>%
  filter(tcell_subset == "Tregs") %>%
  ggplot(aes(x=Day, y=value, fill=therapy_before_cart)) +
  geom_boxplot(alpha=.75, color ="grey20", outlier.shape = NA) +
  geom_point(stat="identity", alpha=.5, color="black", size=1, position = position_jitterdodge(jitter.width = .15), show.legend = F) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45, hjust=1, vjust=1)
  ) +
  facet_grid(.~paste0("Treg")) +
  xlab("Time point") +
  ylab("# cells 10^6/ml") +
  scale_fill_manual(name="Bridging", values=anno_nejm) +
  scale_y_log10(expand = c(0.1,0)) +
  stat_compare_means(label="p.format")

ggsave2(
  "figures/main/figure_5/treg_boxplots.png",
  tregs_bridging,
  width = 103, height = 84, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 5B : Bispecific Ab vs others abundance
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Get proportions
recode_tcells = "'Tregs'='total Tregs'; 'Memory Tregs' = 'memory Tregs'; 'CD8 Thymusemigrants' = 'CD8 thymic emigrants'; 'CD8 Effektorzellen'= 'CD8 effector'; 'CD8 Eff-Memory' = 'CD8 effector memory';
  'CD8 central Memory' = 'CD8 central memory'; 'CD4 Thymusemi' = 'CD4 thymic emigrants'; 'CD4 Effektorzellen'='CD4 effector'; 'CD4 Eff-Memory'='CD4 effector memory'; 'CD4 central Memory'='CD4 central memory'"

t_cell_perc_col <- grep("?%", colnames(tcell.df), value=T)

# Perform Wilcoxon tests between BT groups on all time points
tcell.df_long = tcell.df %>%
  pivot_longer(cols = c(4:ncol(tcell.df)), names_to = "tcell_subset") %>%
  mutate(group = case_when(
    grepl("CD4", tcell_subset) ~ "CD4",
    grepl("CD8", tcell_subset) ~ "CD8",
    grepl("Treg", tcell_subset) ~ "Treg"
  ))

t_cell_df = tcell.df_long %>%
  filter(tcell_subset %in% t_cell_perc_col) %>%
  left_join(clin_data[,c("patient_id", "bridging_bin")]) %>%
  mutate(comp = paste0(tcell_subset, "_", Day), bridging_bin = factor(bridging_bin, levels = c("Other", "Bispecific Ab"))) %>%
  filter(!is.na(value), group != "Treg")

pval_df = t_cell_df %>%
  group_by(tcell_subset, Day) %>%
  wilcox_test(value ~ bridging_bin) %>%
  mutate(p.signif = pval_to_signif_code(p), celltype_short = case_when(
    grepl("CD4", tcell_subset) ~ "CD4",
    grepl("CD8", tcell_subset) ~ "CD8"
  ),
  tcell_subset = car::recode(gsub(" %", "", tcell_subset), recode_tcells))
pval_df$tcell_subset

by_group = t_cell_df %>% 
  group_by(Day, tcell_subset, bridging_bin) %>% 
  summarise(
    count = n(),
    avg_val = mean(value) 
  )

FC.df = by_group %>%
  summarize(FC = avg_val[2] / avg_val[1] ) %>%
  mutate(z_score = scale(FC, center=T, scale=T)) %>%
  ungroup() %>%
  mutate( tcell_subset = car::recode(gsub(" %", "", tcell_subset), recode_tcells))

res = left_join(pval_df, FC.df, by=c("Day", "tcell_subset"))

cd4_diff = res %>%
  filter(celltype_short == "CD4") %>%
  mutate(tcell_subset = factor(tcell_subset, levels = rev(c("CD4 thymic emigrants", "naive CD4", "CD4 central memory", "CD4 effector memory",
                                                            "CD4 effector", "CD8 thymic emigrants", "naive CD8", "CD8 central memory", "CD8 effector memory", "CD8 effector"))),
           Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                        levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
           )) %>%
  ggplot(aes(x=Day, y=tcell_subset, fill = z_score)) +
  geom_tile(alpha=.9) +
  geom_text( aes(label = p.signif), size=5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank() ) +
  theme( axis.title = element_blank(), plot.title = element_text(face="bold")) +
  scale_fill_scico(palette="vik", begin =.1, end =.9, name = "Z-score\n(Bispecific Ab \nvs Other)", midpoint=0) +
  facet_grid(.~paste0("CD4+ subsets"))

cd8_diff = res %>%
  filter(celltype_short == "CD8") %>%
  mutate(tcell_subset = factor(tcell_subset, levels = rev(c("CD4 thymic emigrants", "naive CD4", "CD4 central memory", "CD4 effector memory",
                                                            "CD4 effector", "CD8 thymic emigrants", "naive CD8", "CD8 central memory", "CD8 effector memory", "CD8 effector"))),
         Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                      levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
         )) %>%
  ggplot(aes(x=Day, y=tcell_subset, fill = z_score)) +
  geom_tile(alpha=.9) +
  geom_text(data = subset(res, p<0.05), aes(label = p.signif), size=5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1) ) +
  theme( axis.title = element_blank(), plot.title = element_text(face="bold")) +
  scale_fill_scico(palette="vik", begin =0, end =1, name = "Z-score\n(Bispecific Ab \nvs Other)", , midpoint=0) +
  facet_grid(.~paste0("CD8+ subsets"))

cd4_diff + cd8_diff + plot_layout(nrow=2, guides="collect")

ggsave2(
  "figures/main/figure_5/t_diff_heatmap.png",
  width = 100, height=90, dpi = 400, units = "mm", bg="white",
  scale=1.3
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 5C : Checkpoints on CARs
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

marker_lvl = c("PD1", "TIGIT", "TIM3", "LAG-3", "VISTA")
day_lvl = c("Day 7", "Day 14", "Day 30", "Day 100")

pd1_long = pd1_car_df %>%
  dplyr::select(-patient_id) %>%
  pivot_longer(c(2:9)) %>%
  filter(!is.na(value)) %>%
  separate_wider_delim(name, "_", names = c("Day", "name")) %>%
  rowwise() %>%
  mutate(Day = factor(Day, levels = day_lvl), line = str_split(name, "\\+")[[1]][1], marker = str_split(name, "\\+")[[1]][3])

checkpoints_long = checkpoint_immu_df %>%
  pivot_longer(c(3:34)) %>%
  filter(grepl("CAR\\+CD[4,8]\\+[A-Z]", name), Day %in% day_lvl) %>%
  rowwise() %>%
  mutate(line = str_split(name, "\\+")[[1]][2], marker = str_split(name, "\\+")[[1]][3])

longitudinal_checkpoint_cars = rbind(checkpoints_long, pd1_long) %>%
  left_join(clin_data[,c("sample_id", "bridging_bin")]) %>%
  mutate(marker = factor(marker, levels = marker_lvl)) %>%
  ggplot(aes(x=Day, y=value, fill=bridging_bin)) +
  geom_boxplot(outlier.shape = NA, linewidth=0.3, alpha=.75) +
  geom_point(color = "grey20", alpha = .6, position = position_jitterdodge(jitter.width = .3), size=.8) +
  scale_fill_manual("Bridging", values=c(anno_nejm[1], Other="grey60")) +
  facet_grid(line~marker) +
  mytheme_grid() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), panel.grid.major.x = element_blank(), legend.position = "bottom", panel.spacing = unit(.75, "lines")) +
  ylab("Proportion (%) of CAR+ cells") +
  xlab("Time point after CAR T cell infusion") +
  scale_y_continuous(expand = c(0.1, 0), breaks=seq(0,100,25)) +
  geom_pwc(label="p.signif", hide.ns=T)

ggsave2(
  "figures/main/figure_5/longitudinal_checkpoint_cars.png",
  longitudinal_checkpoint_cars,
  width = 215, height = 109, dpi = 400, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 5D : PD1 violin plots
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# pvalue (ns as "ns")
pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*", "ns"))))
}

longitudinal_pd1_bite_df <- immune_subset_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")),
    bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other") ) %>% 
  filter(
    # therapy_before_cart == "Bispecific Ab",
    grepl("^CD[0-9]{1}_PD1_%", immune_subset),
    immune_subset != "CD3_PD1_%",
    bridging_bin == "Bispecific Ab"
  ) %>%
  mutate(immune_subset = car::recode(immune_subset, "'CD4_PD1_%'='CD4+PD1+'; 'CD8_PD1_%' = 'CD8+PD1+'"))

longitudinal_pd1_bite_cd4_p <- JonckheereTerpstraTest(Day ~ perc, subset(longitudinal_pd1_bite_df, subset = immune_subset == "CD4+PD1+"))
longitudinal_pd1_bite_cd8_p <- JonckheereTerpstraTest(Day ~ perc, subset(longitudinal_pd1_bite_df, subset = immune_subset == "CD8+PD1+"))
bite_jockheere_df <- data.frame(
  immune_subset = c("CD4+PD1+", "CD8+PD1+"),
  p_signif =  c(pval_to_signif_code(longitudinal_pd1_bite_cd4_p$p.value), pval_to_signif_code(longitudinal_pd1_bite_cd8_p$p.value)),
  bridging_bin = "Bispecific Ab"
)

longitudinal_pd1_other_df <- immune_subset_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")),
    bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other") ) %>% 
  filter(
    grepl("^CD[0-9]{1}_PD1_%", immune_subset),
    immune_subset != "CD3_PD1_%",
    bridging_bin == "Other"
  ) %>%
  mutate(immune_subset = car::recode(immune_subset, "'CD4_PD1_%'='CD4+PD1+'; 'CD8_PD1_%' = 'CD8+PD1+'"))

longitudinal_pd1_other_cd4_p <- JonckheereTerpstraTest(Day ~ perc, subset(longitudinal_pd1_other_df, subset = immune_subset == "CD4+PD1+"))
longitudinal_pd1_other_cd8_p <- JonckheereTerpstraTest(Day ~ perc, subset(longitudinal_pd1_other_df, subset = immune_subset == "CD8+PD1+"))
other_jockheere_df <- data.frame(
  immune_subset = c("CD4+PD1+", "CD8+PD1+"),
  p_signif =  c(pval_to_signif_code(longitudinal_pd1_other_cd4_p$p.value), pval_to_signif_code(longitudinal_pd1_other_cd8_p$p.value)),
  bridging_bin = "Other"
)

# Joint and switched
longitudinal_pd1_df = immune_subset_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")),
    bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other") ) %>% 
  filter(
    grepl("^CD[0-9]{1}_PD1_%", immune_subset),
    immune_subset != "CD3_PD1_%",
  ) %>%
  mutate(immune_subset = car::recode(immune_subset, "'CD4_PD1_%'='CD4+PD1+'; 'CD8_PD1_%' = 'CD8+PD1+'"))

jockheere_df = rbind(bite_jockheere_df, other_jockheere_df)

longitudinal_pd1_violins = ggplot(longitudinal_pd1_df, aes(x=Day, y=perc, fill=bridging_bin)) +
  geom_violin(alpha=.75, draw_quantiles = .5, color ="grey20", show.legend = F) +
  geom_point(color = "grey20", alpha = .6, position = position_jitterdodge(jitter.width = .3), size=.8) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank()
  ) +
  xlab("Time point after CAR T cell infusion") +
  ylab("Proportion (%) of non-transduced T cells") +
  scale_fill_manual("Bridging", values=c(anno_nejm[1], Other="grey60")) +
  geom_bracket(
    data = jockheere_df, aes(label = p_signif),
    xmin = "LA", xmax = "Day 100", y.position = 103,
  ) +
  facet_grid(immune_subset~bridging_bin, scales="free_y") +
  scale_y_continuous(expand = c(.12,.07), breaks=seq(0,100,25))

ggsave2(
  filename="figures/main/figure_5/longitudinal_pd1_violins.png",
  longitudinal_pd1_violins,
  width = 130, height = 140, dpi = 400, bg = "white", units = "mm", scale = 1
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 5E : T cell exhaustion signature on single cell data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

lab = as_labeller((c("Apheresis"="LA", "Late"="Day 30", "Very Late"="Day 100", "CD4 T-Cell"="CD4+", "CD8 T-Cell"="CD8+")))

tcell_exh_sig = sc.t %>%
  filter(TIMEPOINT %in% c("Apheresis", "Late", "Very Late"), celltype_short_2 %in% c("CD8 T-Cell", "CD4 T-Cell")) %>%
  group_by(TIMEPOINT, celltype_short_2, celltype, BRIDGING_BIN) %>%
  summarize(mean = mean(exhaustion_UCell_2)) %>%
  ungroup() %>%
  mutate(scaled = scale(mean)) %>%
  ggplot(aes(x=BRIDGING_BIN, y=celltype, fill = scaled)) +
  geom_tile() +
  facet_grid(celltype_short_2~TIMEPOINT, scales="free", labeller = lab) +
  scale_fill_scico(palette = "vik", begin=0, end=.95, midpoint=0, name="Z-score\nexhaustion") +
  theme(axis.title = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust=1))

ggsave2(
  "figures/main/figure_5/exhaustion_heatmap.png",
  tcell_exh_sig,
  height=110, width=133, units="mm", bg="white", dpi=400, scale=1.2
)
