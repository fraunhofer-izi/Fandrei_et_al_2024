.cran_packages = c("tidyverse", "psych", "cowplot", "DescTools", "ggalluvial", 
                   "ggtext")

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
immune_df = bridging.obj$immune_df
t_cells = bridging.obj$t_cells
cytotox_df = bridging.obj$cytotox

sc.pl = read.csv(
  "data/cpc_data.csv"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S1 : PFS analysis BsAb vs other
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

km_bridging_bin = km_plot(
  clin_data,
  "pfs_months",
  "status",
  "bridging_bin",
  c(anno_nejm[1], Other = "grey60"),
  "Bridging",
  c(10, 0.95),
  c(0, 12),
  break.by=2
)

ggsave2(
  "figures/supplement/km_bridging_bin.png",
  km_bridging_bin,
  width = 140, height = 115, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S2: Correlation of CAR-T cell expansion by PCR and flow cytometry
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

corr_df <- immune_df %>%
  filter(!Day %in% c("Leukapheresis", "Day 0")) %>%
  dplyr::select(patient_id, Day, CD3_CAR_abs) %>%
  mutate(CD3_CAR_abs = as.numeric(CD3_CAR_abs) ) %>%
  left_join(pcr_data, by=c("patient_id", "Day")) %>%
  filter(!is.na(sample))

corr_data = corr_df %>%
  group_by(Day, Product) %>%
  filter(!is.na(CD3_CAR_abs), !is.na(no_copies)) %>%
  summarize(
    COR = stats::cor.test(no_copies, CD3_CAR_abs, method="spearman")$estimate,
    pval = stats::cor.test(no_copies, CD3_CAR_abs, method="spearman")$p.value,
    min_x = min(log10(no_copies+1)), max_x = max(log10(no_copies+1))
  ) %>%
  mutate(x=min_x + (max_x-min_x)/2)

corr_scatter = corr_df %>%
  mutate(Day = factor(Day, levels=c("Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  ggplot(aes(x=log10(no_copies+1), y=log10(CD3_CAR_abs+1), color=Product )) +
  geom_smooth(method='lm', formula= y~x, alpha=.75, se=T, fill="grey90") +
  geom_point(size=1.5, alpha=.7) +
  mytheme_grid(13) +
  theme(
    legend.position = "bottom",
    axis.title.y = element_markdown()
  ) +
  scale_color_manual(name="CAR-T Product", values=pal_product) +
  geom_text(data=corr_data, aes(x=x, y=4.9, label=paste0("r=", signif(COR,2), ", p=", signif(pval, 2))), color="black") +
  scale_y_continuous(expand=c(0.12,0)) +
  ylab(bquote("\\# log(x+1) CD3+CAR+ 10<sup>6</sup>/mL")) +
  xlab("# log(x+1) copies per 1000 cells") +
  facet_grid2(Product~Day, scales="free", independent = "all")

ggsave2(
  "figures/supplement/corr_scatter.png",
  corr_scatter,
  width = 150, height = 100, dpi = 400, bg = "white", units = "mm", scale = 1.5
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S3: Impact of bridging therapy on lymphodepletion
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ldp_boxes = t_cells_abs %>%
  mutate(Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"), levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  filter(group != "Treg") %>%
  ggplot(aes(x=Day, y = value, fill = bridging_bin)) +
  geom_boxplot(alpha=.8, outlier.size=.4, linewidth=.4) +
  facet_wrap(~ tcell_subset, scales = "free_y", ncol=2) +
  scale_fill_manual(name="Bridging", values=c(anno_nejm[[1]], "grey60")) +
  mytheme_grid(11) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), legend.position="bottom", panel.grid.major.x = element_blank(),
        axis.title.y = element_markdown()) +
  geom_pwc(label="p.signif", hide.ns="p") +
  scale_y_continuous(expand=c(0.2,0), trans="log10") +
  ylab("Cell count 10<sup>6</sup>/mL") +
  xlab("Time point")

ggsave2(
  "figures/supplement/ldp_boxes.png",
  ldp_boxes,
  width = 85, height = 130, dpi = 400, bg = "white", units = "mm", scale = 1.7
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S4: CAR-T expansion by flow cytometry using absolute counts in the
# bridging therapy groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pl_car_abs = immune_df %>%
  left_join(clin_data[,c("patient_id", "therapy_before_cart")]) %>%
  filter(!Day %in% c("Leukapheresis", "Day 0")) %>%
  ggplot(aes(x=therapy_before_cart, y=as.numeric(`CD3_CAR_abs`), fill=therapy_before_cart)) +
  geom_boxplot(alpha=.8, outlier.shape = NA) +
  mytheme_grid(13) +
  geom_point(color="grey10", position = position_jitter(width=0.3), alpha=.7, size=1.2) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown()
  ) +
  facet_grid(.~Day, scales="free_y") +
  ylab("Number of CAR-T cells 10<sup>6</sup>/mL") +
  scale_fill_manual(name = "Bridging", values=anno_nejm) +
  stat_compare_means(label="p.format",hide.ns=F, size=3.5, hjust=0.3, vjust=0) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, label.size=5.5) +
  scale_y_log10()

ggsave2(
  "figures/supplement/total_car_expansion.png",
  pl_car_abs,
  width = 120, height = 65, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S5: Proportion of CAR-T cells in samples extracted on day 7 after CAR-T
# infusion for in vitro cytotoxicity assay
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cytotox_prop_cart =  cytotox_df %>%
  filter(name == "percentage_cart") %>%
  ggplot(aes(x=therapy_before_cart, y=value, fill=therapy_before_cart)) +
  geom_boxplot(alpha=.75, outlier.shape = NA) +
  geom_point(color="grey10", position = position_jitterdodge(jitter.width=0.4), alpha=.7, size=1.3) +
  mytheme_grid(13) +
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face="bold")
  ) +
  scale_fill_manual("Bridging", values=anno_nejm) +
  stat_compare_means(label = "p.format", vjust=1) +
  geom_pwc(label="p.signif", hide.ns="p") +
  ylab("Proportion (%) of CAR-T cells") +
  ggtitle("Proportion of CAR-T cells in expanded PB samples collected\n7d after infusion")

ggsave2(
  "figures/supplement/cytotox_prop_cart.png",
  cytotox_prop_cart,
  height=75, width=100, units = "mm", dpi=400, scale=1.6, bg="white"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S5: Impact of bridging therapy on Treg subsets
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

treg_subset_violins <- treg_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
    )) %>% 
  filter(tcell_subset != "Tregs%") %>%
  left_join(clin_data[,c("patient_id", "therapy_before_cart")]) %>%
  ggplot(aes(x=therapy_before_cart, y=perc*100, fill=therapy_before_cart)) +
  geom_boxplot(alpha=.75, outlier.shape = NA, color ="grey20") +
  geom_point(stat="identity", alpha=.5, color="black", size=1, position = position_jitterdodge(jitter.width = .4), show.legend = F) +
  mytheme_grid(12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45, hjust=1),
    axis.title.x = element_blank()
  ) +
  facet_grid2(tcell_subset~Day, scales = "free_y") +
  ylab("%  of total Treg") +
  scale_fill_manual(name="Bridging", values=anno_nejm) +
  scale_y_continuous(expand = c(0.05, 0), breaks=pretty_breaks(n=4)) +
  stat_compare_means(label="p.format", vjust=-3) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, vjust=0.5, label.size=5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

ggsave2(
  "figures/supplement/treg_subsets.png",
  treg_subset_violins,
  width=185, height=120, dpi=400, units="mm", bg="white", scale=1.3
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S6: T cell subsets according to remission status before CAR-T
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Prepare
recode_tcells <- "'CD4 Thymusemi %' = 'CD4+ thymic emigrants'; 'naive CD4 %' = 'Naive CD4+'; 'CD4 Eff-Memory %' = 'CD4+ effector memory';
    'CD4 central Memory %' = 'CD4+ central memory'; 'CD4 Effektorzellen %' = 'CD4+ effector'; 'CD8 Thymusemigrants %' = 'CD8+ thymic emigrants';
    'naive CD8 %' = 'Naive CD8+'; 'CD8 Eff-Memory %' ='CD8+ effector memory'; 'CD8 central Memory %' = 'CD8+ central memory';
    'CD8 Effektorzellen %' = 'CD8+ effector'; 'Tregs %' = 'Tregs';  'naive Tregs %' = 'Naive Tregs';
    'Memory Tregs %'= 'Memory Tregs'"

tcell.df_long = tcell.df %>%
  pivot_longer(cols = c(4:ncol(tcell.df)), names_to = "tcell_subset") %>%
  mutate(group = case_when(
    grepl("CD4", tcell_subset) ~ "CD4",
    grepl("CD8", tcell_subset) ~ "CD8",
    grepl("Treg", tcell_subset) ~ "Treg"
  ))

pretty_unexpanded <- function(x,n=3){
  #expand_limit
  if(x[1]<=0){
    r2 <- x + c(-x[1],x[1])
  }else{
    r2 <- x + c((x[2]-x[1])*0.04545455,-(x[2]-x[1])*0.04545455)
  }
  pout <-  pretty(r2,n) 
  pout
}

pl_status_before_cart_cd4 = tcell.df_long %>%
  filter(group == "CD4", grepl("?%", tcell_subset)) %>%
  mutate(tcell_subset = car::recode(tcell_subset, recode_tcells),
         Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
           levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  left_join(clin_data[,c("patient_id", "status_before_cart")]) %>%
  ggplot(aes(x=status_before_cart, y=100*value, fill=status_before_cart)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size=.8, color="grey20", alpha=.8) +
  facet_grid2(tcell_subset ~ Day, scales="free_y") +
  mytheme_grid(12) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        plot.tag = element_text(face="bold", family = "serif", size = rel(2))) +
  xlab("Remission status before CAR-T") +
  scale_fill_manual(name="", values=palette_response) +
  geom_pwc(label="p.signif", hide.ns="p") +
  stat_compare_means(label="p.format", vjust=-2.25, ) +
  scale_y_continuous(expand = c(0.36, 0), limits = c(0, NA), breaks = pretty_unexpanded) +
  ylab("Proportion (%) of CD4+ T cells") +
  labs(tag = 'A)')

ggsave2(
  "figures/supplement/pl_status_before_cart_cd4.png",
  pl_status_before_cart_cd4,
  height=205, width=160, units="mm", bg="white", dpi=400, scale=1.4
)

pl_status_before_cart_cd8 = tcell.df_long %>%
  filter(group == "CD8", grepl("?%", tcell_subset)) %>%
  mutate(tcell_subset = car::recode(tcell_subset, recode_tcells),
         Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                      levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100"))) %>%
  left_join(clin_data[,c("patient_id", "status_before_cart")]) %>%
  ggplot(aes(x=status_before_cart, y=100*value, fill=status_before_cart)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size=.8, color="grey20", alpha=.8) +
  facet_grid2(tcell_subset ~ Day, scales="free_y") +
  mytheme_grid(12) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        plot.tag = element_text(face="bold", family = "serif", size = rel(2))) +
  xlab("Remission status before CAR-T") +
  scale_fill_manual(name="", values=palette_response) +
  geom_pwc(label="p.signif", hide.ns="p") +
  stat_compare_means(label="p.format", vjust=-2.25) +
  scale_y_continuous(expand=c(0.36, 0), breaks=pretty_unexpanded) +
  ylab("Proportion (%) of CD8+ T cells") +
  labs(tag = 'B)')

ggsave2(
  "figures/supplement/pl_status_before_cart_cd8.png",
  pl_status_before_cart_cd8,
  height=205, width=160, units="mm", bg="white", dpi=400, scale=1.4
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S7: T cell subset alluvials
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# T cell subsets
recode_tcells <- "'CD4 Thymusemi %' = 'CD4+ thymic emigrants'; 'naive CD4 %' = 'Naive CD4+'; 'CD4 Eff-Memory %' = 'CD4+ effector memory';
    'CD4 central Memory %' = 'CD4+ central memory'; 'CD4 Effektorzellen %' = 'CD4+ effector'; 'CD8 Thymusemigrants %' = 'CD8+ thymic emigrants';
    'naive CD8 %' = 'Naive CD8+'; 'CD8 Eff-Memory %' ='CD8+ effector memory'; 'CD8 central Memory %' = 'CD8+ central memory';
    'CD8 Effektorzellen %' = 'CD8+ effector'; 'Tregs %' = 'Tregs';  'naive Tregs %' = 'Naive Tregs';
    'Memory Tregs %'= 'Memory Tregs'"

t_levels = c("CD4+ central memory", "CD4+ effector", "CD4+ effector memory", "CD4+ thymic emigrants",
             "Naive CD4+", "CD8+ central memory", "CD8+ effector", "CD8+ effector memory",
             "CD8+ thymic emigrants", "Naive CD8+")

alluvial_t_cells = t_cells %>%
  pivot_longer(c(4:ncol(t_cells)), names_to = "tcell_subset") %>%
  mutate(group = case_when(
    grepl("CD4", tcell_subset) ~ "CD4",
    grepl("CD8", tcell_subset) ~ "CD8",
    grepl("Treg", tcell_subset) ~ "Treg"
  ),
  Day = factor(ifelse(Day == "Leukapheresis", "LA", gsub("Day ", "", Day)),
               levels=c("LA", "0", "7", "14", "30", "100"))
  ) %>%
  filter(grepl("?%", tcell_subset), group != "Treg") %>%
  mutate(tcell_subset = factor(car::recode(tcell_subset, recode_tcells), levels = t_levels)) %>%
  left_join(clin_data[,c("patient_id", "therapy_before_cart")]) %>%
  group_by(therapy_before_cart, Day, group, tcell_subset) %>%
  summarize(mean = mean(100*value)) %>%
  ggplot(aes(x = Day, y = mean, alluvium = tcell_subset)) +
  geom_alluvium(aes(fill = tcell_subset, colour = tcell_subset), alpha = .9, decreasing = F, color = "grey20") +
  theme( axis.title.x = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
  scale_fill_manual(name = "T cell subset", values = t_pal) +
  facet_grid(group ~ therapy_before_cart, scales = "fixed") +
  ylab("Proportion (%)") +
  guides(fill = guide_legend(ncol=2), color = guide_legend(ncol=2))

ggsave2(
  "figures/supplement/alluvial_t_cells.png",
  alluvial_t_cells,
  width = 115, height = 95, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S8: PD-1 violins - all comparison
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Compare each timepoint
pd1_violins <- immune_subset_df %>%
  mutate( 
    Day = factor(car::recode(Day, "'Leukapheresis' = 'LA'"),
                 levels=c("LA", "Day 0", "Day 7", "Day 14", "Day 30", "Day 100")
    )) %>% 
  filter(
    # therapy_before_cart == "Bispecific Ab",
    grepl("^CD[0-9]{1}_PD1_%", immune_subset),
    immune_subset != "CD3_PD1_%"
  ) %>%
  mutate(immune_subset = car::recode(immune_subset, "'CD4_PD1_%'='CD4+'; 'CD8_PD1_%' = 'CD8+'")) %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill=therapy_before_cart)) +
  geom_violin(alpha=.75, draw_quantiles = .5, color ="grey20") +
  geom_point(color="black", alpha=.5, size=1, position = position_jitterdodge(jitter.width=.35)) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank()
  ) +
  facet_grid(immune_subset~Day, scales="free_y") +
  xlab("Days since intervention") +
  ylab("Proportion (%) of PD-1 expressing cells") +
  scale_fill_manual(name="Bridging", values=anno_nejm) +
  scale_y_continuous(breaks=c(0,25,50,75,100), lim=c(0,105), expand=c(0.2, 0)) +
  stat_compare_means(label="p.format", vjust=-1) +
  geom_pwc(method="wilcox.test", label = "p.signif", hide.ns=T, vjust=0.5, label.size=5) +
  guides(fill = guide_legend(nrow = 2))

ggsave2(
  "figures/supplement/pd1_violins.png",
  pd1_violins,
  width=130, height=90, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. S11: Numbers of circulating plasma cells
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

