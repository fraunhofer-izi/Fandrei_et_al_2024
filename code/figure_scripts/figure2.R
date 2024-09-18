.cran_packages = c("tidyverse", "psych", "cowplot", "ggpubr", "ggh4x")

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
crs_df = bridging.obj$crs

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2A : risk categories
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# labels
freq_df = clin_data %>% 
  group_by(therapy_before_cart) %>%
  summarize(Freq = n())
freq_diag <- data.frame(freq_df %>% mutate(labels = paste0(as.character(therapy_before_cart), "\n", "n= ", Freq)))
freq_list <- freq_diag$labels
freq_list

refractoriness_plot <- clin_data %>%
  group_by(therapy_before_cart, refractoriness_group) %>%
  mutate(
    therapy_before_cart = factor(therapy_before_cart, levels=levels(therapy_before_cart)),
    refractoriness_group = factor(refractoriness_group, levels = c("TCexposed", "TCRRMM", "PentaRRMM"))
  ) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n)) %>%
  ggplot(aes(x=fct_rev(therapy_before_cart), y=perc, fill = refractoriness_group)) +
  geom_bar(stat="identity", alpha=.9, color="black", position = position_stack(reverse=T)) +
  theme_light() + gtheme(13) + 
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    strip.background =element_rect(fill="grey45"),
    strip.text = element_text(colour = 'white'),
    axis.line=element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_fill_manual( values=refractoriness_col) +
  facet_grid(.~paste0("Refractoriness")) +
  xlab("Bridging") + ylab("") +
  scale_x_discrete(labels = rev(freq_list)) +
  coord_flip() +
  guides(fill=guide_legend(nrow=3, title = NULL, reverse = TRUE))

fisher.test(clin_data$bridging_bin, clin_data$refractoriness_group)

# Number of high-risk cytogenetics
cytogen_df <- clin_data %>%
  dplyr::select(c(patient_id, therapy_before_cart, `t(4;14)`, del_17p, gain_1q)) %>% 
  mutate(across(c(3:5), function(x) as.numeric(paste(x)))) %>%
  mutate(sum = `t(4;14)` + del_17p + gain_1q) %>%
  mutate(high_risk = ifelse(sum>0, 1, 0)) %>%
  left_join(clin_data[,c("patient_id", "bridging_bin")], by="patient_id")

table(cytogen_df$therapy_before_cart, cytogen_df$sum)
chisq.test(cytogen_df$therapy_before_cart, cytogen_df$sum)

table(cytogen_df$bridging_bin, cytogen_df$high_risk)
fisher.test(cytogen_df$bridging_bin, cytogen_df$high_risk)

cytogen_plot <- cytogen_df %>% 
  group_by(therapy_before_cart, sum) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n), sum = factor(sum, levels = c(0,1,2,3))) %>%
  ggplot(aes(x=fct_rev(therapy_before_cart), y=perc, fill = sum)) +
  geom_bar(position = position_stack(reverse=T), stat="identity", alpha=.9, color="black") +
  theme_light() + gtheme(13) + 
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    strip.background =element_rect(fill="grey45"),
    strip.text = element_text(colour = 'white'),
    axis.line=element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_x_discrete(labels = rev(freq_list)) +
  scale_fill_manual( values=col.refr) +
  facet_grid(.~paste0("High-risk Cytogen")) +
  xlab("Bridging") + ylab("") +
  coord_flip() +
  guides(fill=guide_legend(nrow=4, title = "Number of \nHR Cytogen", title.position = "left", reverse = TRUE))

# MyCARe
refractoriness_df = bridging.obj$refractoriness %>% 
  rename(len_refr  = "Lenalidomide") %>%
  mutate(len_refr = ifelse(len_refr == "refractory", 1, 0))

mycare_df <- clin_data %>%
  dplyr::select(patient_id, therapy_before_cart, EMD, `Ferritin>norm`) %>%
  rename(fer = "Ferritin>norm") %>%
  left_join(cytogen_df[,c("patient_id", "high_risk")]) %>%
  left_join(refractoriness_df[,c("patient_id", "len_refr")]) %>% 
  mutate(fer = as.numeric(paste(fer)), EMD = as.numeric(paste(EMD))) %>%
  mutate(score = rowSums(.[3:6])) %>% 
  mutate(mycare = case_when(
    score %in% c(0, 1) ~ "Low",
    score %in% c(2, 3) ~ "Intermediate",
    score == 4 ~ "High"
  )) %>% 
  left_join(clin_data[,c("patient_id", "bridging_bin")], by="patient_id")

fisher.test(mycare_df$therapy_before_cart, mycare_df$mycare)

mycare_df %>%
  group_by(bridging_bin, mycare) %>%
  summarize(n = n()) %>%
  mutate(prop = n/sum(n))

fisher.test(mycare_df$bridging_bin, mycare_df$mycare)

mycare_plot <- mycare_df %>%
  group_by(therapy_before_cart, mycare) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n), mycare = factor(mycare, levels = c("Low", "Intermediate", "High"))) %>%
  ggplot(aes(x=fct_rev(therapy_before_cart), y=perc, fill = mycare)) +
  geom_bar(position = position_stack(reverse=T), stat="identity", alpha=.9, color="black") +
  theme_light() + gtheme(13) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    strip.background =element_rect(fill="grey45"),
    strip.text = element_text(colour = 'white'),
    axis.line=element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_fill_manual( values=col.refr[c(1,2,3)] ) +
  facet_grid(.~paste0("MyCARe")) +
  xlab("Bridging") + ylab("") +
  coord_flip() +
  guides(fill=guide_legend(nrow=3, title = "MyCARe \nScore", title.position = "left", reverse = TRUE))

riss_plot <- clin_data %>%
  filter(!is.na(`R-ISS`)) %>%
  group_by(therapy_before_cart, `R-ISS`) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n), `R-ISS` = factor(`R-ISS`, levels = c(1, 2, 3))) %>%
  ggplot(aes(x=fct_rev(therapy_before_cart), y=perc, fill = `R-ISS`)) +
  geom_bar(position = position_stack(reverse=T), stat="identity", alpha=.9, color="black") +
  theme_light() + gtheme(13) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    strip.background =element_rect(fill="grey45"),
    strip.text = element_text(colour = 'white'),
    axis.line=element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_fill_manual( values=col.refr[c(1,2,3)] ) +
  facet_grid(.~paste0("R-ISS")) +
  xlab("Bridging") + ylab("") +
  coord_flip() +
  guides(fill=guide_legend(nrow=3, title = "R-ISS \nScore", title.position = "left", reverse = TRUE))

risk_panel = refractoriness_plot + cytogen_plot + mycare_plot + riss_plot +
  plot_layout(nrow = 1) +
  plot_annotation(title = "Risk categories in patients with different bridging therapies") &
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

risk_panel_f = gridExtra::grid.arrange(patchworkGrob(risk_panel), bottom = "Proportion of patients")

ggsave2(
  "figures/main/figure_2/risk_panel.png",
  plot = risk_panel_f,
  width = 210, height = 80, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2B : ORR of bridging therapy
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

clin_data %>%
  group_by(all_groups_before_cart) %>%
  summarize(n=n()) %>% 
  mutate(perc=n/sum(n)*100)

clin_data %>%
  group_by(bridging_bin, status_before_cart) %>%
  summarize(n=n()) %>% 
  mutate(perc=n/sum(n)*100)

orr_box <- clin_data %>%
  group_by(all_groups_before_cart) %>%
  summarize(n=n()) %>% 
  mutate(perc=n/sum(n)*100) %>%
  filter(!all_groups_before_cart %in% c("SD", "PD")) %>%
  mutate(all_groups_before_cart = factor(all_groups_before_cart, levels=rev(c("CR", "VGPR", "PR"))) ) %>%
  ggplot(aes(y=perc, x="", fill=all_groups_before_cart)) +
  geom_bar(position="stack", stat='identity', color="black") +
  geom_text(aes(label=paste0(all_groups_before_cart, " (", round(perc, 0), " %)")), position = position_stack(vjust = 0.5),
            color="white", size=5, fontface='bold') +
  coord_trans(ylim=c(0,70)) +
  theme_classic() + gtheme(17) + 
  ylab("Response to BT in %") +
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),axis.text.x=element_blank(), axis.title.x = element_blank()) +
  # scale_fill_manual("Response groups", values=c(derain[4:6])) +
  scale_fill_manual("Response groups", values=rev(c("grey30", "grey50", "grey70"))) +
  guides(fill="none") +
  annotate("text", y=65, x=1, label="ORR: 56%", size=5) +
  scale_y_continuous(position = "left")

ggsave2(
  "figures/main/figure_2/orr_box.png",
  plot = orr_box,
  width = 50, height = 75, dpi = 400, bg = "white", units = "mm", scale = 1.4
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2C : Survival - Remission before CAR-T
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

km_response_to_bridging = km_plot(
  clin_data,
  "pfs_months",
  "status",
  "status_before_cart",
  palette_response,
  "Response to BT",
  c(9, .9),
  limits=c(0,12),
  conf = F,
  break.by=3
)

ggsave2(
  "figures/main/figure_2/km_response_to_bridging.png",
  km_response_to_bridging,
  width = 130, height = 130, dpi = 400, bg = "white", units = "mm", scale = 1.3
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2D: Response after bridging bar plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Binary
clin_data <- clin_data %>%
  mutate(
    bridging_binary = ifelse(therapy_before_cart == "Bispecific Ab", "BITE", "Other"),
    response_bin_before = ifelse(status_before_cart == "SD/PD", "no", "yes")
  )

remission_before_cart_df <- clin_data %>%
  mutate(therapy_before_cart = factor(therapy_before_cart, levels=rev(c("Bispecific Ab", "CD38", "Chemotherapy", "SLAMF7")))) %>%
  group_by(therapy_before_cart, status_before_cart) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n)*100, status = factor(status_before_cart, levels=rev(levels(status_before_cart))))

remission_before_cart_plot <- remission_before_cart_df %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill = fct_rev(status_before_cart))) +
  geom_bar(position="stack", stat="identity", alpha=.75, color="black") +
  theme_classic() + gtheme(18) + theme(axis.line=element_blank(), axis.ticks.y = element_blank()) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_text(size = 14)
  ) +
  scale_fill_manual("Remission before CAR-T", values=palette_response) +
  xlab("Bridging therapy") + ylab("%") +
  guides(fill = guide_legend(reverse = T, nrow = 1, title = NULL)) +
  scale_x_discrete(labels = rev(freq_list)) +
  coord_flip()

orr_data = clin_data %>%
  mutate(therapy_before_cart = factor(therapy_before_cart, levels=rev(c("Bispecific Ab", "CD38", "Chemotherapy", "SLAMF7")))) %>%
  mutate(overall_response = ifelse(status_before_cart=="SD/PD", 0, "ORR")) %>%
  group_by(therapy_before_cart, overall_response) %>%
  summarize(n=n()) %>%
  mutate(orr = n/sum(n)*100) %>%
  filter(overall_response != 0)

clin_data = clin_data %>%
  mutate(overall_response = ifelse(status_before_cart=="SD/PD", 0, "ORR"))

p.orr = fisher.test(table(clin_data$bridging_bin, clin_data$overall_response))

orr_plot = ggplot(orr_data, aes(x=therapy_before_cart, y=orr, fill=overall_response)) +
  geom_bar(position="stack", stat="identity", alpha=.75, color="black", show.legend = F) +
  theme_classic() + gtheme(18) + theme(axis.line=element_blank(), axis.ticks.y = element_blank()) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1)
  ) +
  scale_fill_manual("", values="grey20") +
  xlab("Bridging therapy") + ylab("% ORR") +
  geom_text(inherit.aes = F, aes(label=paste0("P=", signif(p.orr$p.value, 1)), x=1, y=70), size=4) +
  coord_flip() 

response_bridging_g <- remission_before_cart_plot + orr_plot + plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "top")

ggsave2(
  "figures/main/figure_2/response_bridging.png",
  plot = response_bridging_g,
  width = 100, height = 65, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Duration of response
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

dur_df = clin_data %>% filter(`30 day response` != "SD/PD")
survfit(Surv(pfs_months, status) ~ 1, dur_df)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2D : Inflammation markers violins
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

crs_long <- crs_df[,c(1, 5:57)] %>% 
  pivot_longer(!sample_id, names_to = "crs_marker", values_to = "value") %>%
  left_join(crs_df[,c("sample_id", "crs_group")], by="sample_id") %>%
  left_join(clin_data[,c("sample_id", "Product", "therapy_before_cart")], by="sample_id") %>%
  separate(crs_marker, c("marker", "Day"), sep="_", remove=F)

# Time of measurements
crs_long %>%
  filter(!is.na(value)) %>%
  group_by(sample_id, marker) %>%
  slice_tail(n = 1) %>%
  filter(marker == "CRP") %>%
  mutate(day = as.numeric(gsub("d", "", Day))) %>%
  ungroup() %>% 
  summarize(median = median(day, na.rm=F), min = min(day), max=max(day))

crs_long %>%
  filter(!is.na(value)) %>%
  group_by(sample_id, marker) %>%
  slice_tail(n = 1) %>%
  filter(marker == "PCT") %>%
  mutate(day = as.numeric(gsub("d", "", Day))) %>%
  ungroup() %>% 
  summarize(median = median(day, na.rm=F), min = min(day), max=max(day))

crs_long %>%
  filter(!is.na(value)) %>%
  group_by(sample_id, marker) %>%
  slice_tail(n = 1) %>%
  filter(marker == "il6") %>%
  mutate(day = as.numeric(gsub("\\-.+", "", gsub("d", "", Day)))) %>%
  ungroup() %>%
  summarize(median = median(day, na.rm=F), min = min(day), max=max(day))

inflammation_violins <- crs_long %>% 
  mutate(marker = car::recode(marker, "'il6'='IL6 (pg/mL)'; 'CRP'='CRP (mg/L)'; 'PCT'='PCT (mg/L)'", as.factor=T, levels=c("CRP (mg/L)", "PCT (mg/L)", "IL6 (pg/mL)"))) %>%
  filter(!is.na(value)) %>%
  group_by(marker, sample_id) %>%
  summarize(max_crp = max(value)) %>%
  left_join(crs_df[,c("sample_id", "therapy_before_cart")]) %>%
  ggplot(aes(x=therapy_before_cart, y=max_crp, fill=therapy_before_cart)) + 
  geom_violin(alpha=.75, draw_quantiles = .5) +
  geom_jitter(width=0.3, alpha=.5, show.legend = F) +
  mytheme_grid(13) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_fill_manual(name="Bridging", values=anno_nejm) +
  facet_wrap(.~marker, scales="free", nrow=1) +
  stat_compare_means(label="p", hjust=0.2, vjust=-1, size=5) +
  geom_pwc(label="p.signif", hide.ns="p") +
  scale_y_log10(expand=c(0.15,0)) +
  ylab("Maximal level post-infusion") +
  guides(fill = guide_legend(nrow=2))

ggsave2(
  "figures/main/figure_2/inflammation_markers_violins.png",
  plot = inflammation_violins,
  width = 115, height = 75, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2E: : CRS groups
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

crs_groups_bar_1 <- crs_df %>%
  mutate(
    therapy_before_cart=factor(therapy_before_cart, levels=rev(c("Bispecific Ab", "Chemotherapy", "CD38", "SLAMF7"))),
    crs_group = factor(
      car::recode(crs_group, "'no'='no CRS'; 'yes-Toci'='CRS-Toci'; 'yes+Toci'='CRS+Toci'"),
      levels=c('no CRS', 'CRS-Toci', 'CRS+Toci')
  )) %>%
  group_by(therapy_before_cart, crs_group) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n)*100) %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill = crs_group)) +
  geom_bar(position = position_stack(reverse=T), stat="identity", alpha=.8, color="black") +
  theme_classic() + gtheme(18) + theme(axis.line=element_blank(), axis.ticks.y = element_blank()) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    axis.text.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "top"
  ) +
  scale_fill_manual(name="CRS group", values=col.refr[c(1,2,3)]) +
  xlab("Bridging") + ylab("%") +
  guides(fill = guide_legend(title.position = "left", nrow = 2, reverse = T)) +
  scale_x_discrete(labels = rev(freq_list)) +
  coord_flip()

crs_stats <- crs_df %>%
  mutate(
    crs_group = factor(crs, levels=rev(c("no", "yes"))),
    bridging_bin = ifelse(therapy_before_cart == "Bispecific Ab", "Bispecific Ab", "Other")
  )
fisher.test(crs_stats$bridging_bin, crs_stats$crs_group)

crs_groups_bar_2 <- crs_df %>%
  mutate(
    therapy_before_cart=factor(therapy_before_cart, levels=rev(c("Bispecific Ab", "Chemotherapy", "CD38", "SLAMF7")))) %>%
  group_by(therapy_before_cart, crs) %>%
  mutate(crs_group = factor(crs, levels=rev(c("no", "yes")))) %>%
  summarize(n = n()) %>%
  mutate(perc = n/sum(n)*100) %>%
  filter(crs=="yes") %>%
  ggplot(aes(x=therapy_before_cart, y=perc, fill = crs)) +
  geom_bar(position = "stack", stat="identity", alpha=.8, color="black", show.legend = F) +
  theme_classic() + gtheme(18) + theme(axis.line=element_blank(), axis.ticks.y = element_blank()) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "top"
  ) +
  ylab("% CRS") +
  scale_fill_manual("",values="grey20") +
  coord_flip() +
  annotate("label",label="P=0.1",x=4,y=80,size=5,label.size = NA)

crs_plot <- crs_groups_bar_1 + crs_groups_bar_2 + plot_layout(widths = c(2, 1), guides = "collect") &
  theme(legend.position = "top")

crs_df %>%
  left_join(clin_data[,c("patient_id", "bridging_bin")]) %>%
  group_by(bridging_bin, crs_group) %>%
  summarize(n = n()) %>%
  mutate(prop=n/sum(n))

ggsave2(
  "figures/main/figure_2/crs_plot.png",
  plot = crs_plot,
  width = 100, height = 68, dpi = 400, bg = "white", units = "mm", scale = 1.6
)
