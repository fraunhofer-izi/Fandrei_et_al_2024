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

# Single cell data
sc.t = read.csv(
  "data/tcell_data.csv"
)

sc.pl = read.csv(
  "data/cpc_data.csv"
)

div_bridging = read.csv(
  "data/div_bridging.csv"
)

occ.rep = read.csv(
  "data/occ.rep.csv"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6A : Shannon diversity
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

pl_shannon = div_bridging %>%
  mutate(Day = gsub("PBMC \\| ", "", GROUP)) %>%
  mutate(Day = factor(ifelse(Day == "Apheresis", "LA", Day), levels=c("LA", "Day 30", "Day 100")) ) %>%
  ggplot(aes(x=BRIDGING_BIN, y=shannon, fill=BRIDGING_BIN)) +
  geom_boxplot(outlier.shape=NA, fatten=1.5, linewidth=.2) +
  stat_compare_means(label = "p.format", label.x =.8, label.y.npc=.9, size=3.5) +
  scale_fill_manual(values=c(`Bispecific Ab`=anno_nejm[[1]], Other = "grey60")) +
  geom_point(aes(group = BRIDGING_BIN), size=.9, color="grey30", alpha=.7, position=position_jitterdodge()) +
  facet_grid(LIN ~ Day) +
  ylab("Shannon diversity") +
  mytheme_grid(13) +
  theme(
    panel.spacing = unit(.75, "lines"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    axis.title.x = element_blank()
  )

ggsave2(
  "figures/main/figure_6/shannon_bsab_vs_others.png",
  pl_shannon,
  width=90, height=93, dpi=400, bg="white", scale=1.4, units="mm"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6B : Proportions of singleton T cell clones
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

lab = as_labeller(c(Apheresis = "LA", Late = "Day 30", `Very Late` ="Day 100", CD4 = "CD4", CD8="CD8"))

pL_singletons = sc.t %>%
  filter(TIMEPOINT %in% c("Apheresis", "Late", "Very Late")) %>%
  group_by(PATIENT_ID, BRIDGING_BIN, TIMEPOINT, celltype_short_3, cloneSize) %>%
  summarize(n = n()) %>%
  mutate(prop = 100*n/sum(n), celltype_short_3 = gsub(" T-Cell", "", celltype_short_3)) %>%
  filter(cloneSize == "Single (0 < X <= 1)", celltype_short_3 %in% c("CD4", "CD8")) %>%
  ggplot(aes(x=BRIDGING_BIN, y=prop, fill=BRIDGING_BIN)) +
  geom_boxplot(outlier.shape = NA, fatten = 1.5, linewidth=.2) +
  stat_compare_means(label = "p.format", label.x =.8, label.y.npc=.9, size=3.5) +
  scale_fill_manual(values=c(`Bispecific Ab`=anno_nejm[[1]], Other = "grey60")) +
  geom_point(aes(group = BRIDGING_BIN), size=.9, color="grey30", alpha=.7, position=position_jitterdodge()) +
  facet_grid(celltype_short_3 ~ TIMEPOINT, labeller = lab) +
  ylab("% Singletons") +
  mytheme_grid(13) +
  theme(
    panel.spacing = unit(.75, "lines"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    axis.title.x = element_blank()
  )

ggsave2(
  "figures/main/figure_6/singletons_bsab_vs_others.png",
  pL_singletons,
  width=90, height=93, dpi=400, bg="white", scale=1.4, units="mm"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6C : Percent T cell subtypes in clonotype groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

occ.rep = occ.rep %>% mutate(Day = factor(ifelse(TIME == "Apheresis", "LA", TIME), levels=c("LA", "Day 30", "Day 100")) )

pl_occ_rep = ggplot(occ.rep, aes(x=GROUP, y=value, fill=cloneSize)) +
  geom_col(position="fill") +
  facet_grid(Day ~ celltype) +
  scale_fill_manual(values=clono.col) +
  scale_y_continuous(breaks=c(0, 0.5,1)) +
  mytheme(10) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size=rel(1))
  ) +
  ylab("Proportion of cells") +
  guides(fill=guide_legend(nrow=2, title = "Clone size"))

ggsave2(
  "figures/main/figure_6/pl_occ_rep.png",
  pl_occ_rep,
  height=100, width=180, dpi=400, bg="white", scale=1.3, units="mm"
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 6E : BCMA expression
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sc.pl = sc.pl %>%
  mutate(
    Day = factor(Day, levels = c("Pre Bite", "LA", "Day 30", "Day 100")),
    BRIDGING_BIN = factor(BRIDGING_BIN, levels = c("Teclistamab", "Other"))
  )

lab = as_labeller(c("Pre Bite" = "pre BsAb", "LA" = "LA", "Day 30" = "Day 30", "Day 100" = "Day 100"))

cpc_rna = sc.pl %>%
  filter(TYPE=="TNFRSF17_RNA") %>%
  ggplot(aes(x=BRIDGING_BIN, y=EXPRS, fill=BRIDGING_BIN)) +
  geom_violin(linewidth=.2, scale="width", width=.5, alpha=.6) +
  geom_jitter(shape=".", alpha=.6, width=.15, show.legend=F, color="#555555") +
  stat_summary(fun="median", geom="crossbar", width=.4, show.legend = F, lwd=.35) +
  scale_fill_manual(name="Bridging", values=c("#994455", "grey60")) +
  mytheme_grid(14) +
  facet_wrap(~Day, nrow=1, labeller = lab) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), plot.title = element_text(face="bold.italic"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  ylab("RNA Expression") +
  ggtitle("TNFRSF17") +
  scale_y_continuous(expand=c(.1,0)) +
  stat_compare_means(label="p.format", size=4, vjust=.25)

cpc_adt = sc.pl %>%
  filter(TYPE=="TNFRSF17_ADT") %>%
  ggplot(aes(x=BRIDGING_BIN, y=EXPRS, fill=BRIDGING_BIN)) +
  geom_violin(linewidth=.2, scale="width", width=.5, alpha=.6) +
  geom_jitter(shape=".", alpha=.6, width=.15, show.legend=F, color="#555555") +
  stat_summary(fun="median", geom="crossbar", width=.4, show.legend = F, lwd=.35) +
  scale_fill_manual(name="Bridging", values=c("#994455", "grey60")) +
  mytheme_grid(14) +
  facet_wrap(~Day, nrow=1, labeller = lab) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), plot.title = element_text(face="bold"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  ylab("ADT Expression") +
  ggtitle("BCMA") +
  scale_y_continuous(expand=c(.1,0)) +
  stat_compare_means(label="p.format", size=4, vjust=.25)
  
cpc_violins = cpc_rna + cpc_adt + plot_layout(guides="collect", nrow=1) &
  theme(legend.position = "bottom")

ggsave(
  "figures/main/figure_6/cpc_violins.png",
  cpc_violins,
  height=62, width=180, bg="white", dpi=400, units="mm", scale=1.6
)  
  