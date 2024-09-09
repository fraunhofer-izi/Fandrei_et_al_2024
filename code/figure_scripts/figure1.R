.cran_packages = c("tidyverse", "scales", "psych", "cowplot", "ggsci",
                   "MetBrewer", "ggstar", "ggpubr", "ggh4x", "ggnewscale",
                    "rstatix", "swimplot", "patchwork")


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

clin_data =  bridging.obj$clin_data
bite_admin = bridging.obj$bite_admin
apheresis_oos_df = bridging.obj$bite_leukapheresis

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Follow-up
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

surv_formula = as.formula(paste0("Surv(os_months, death)~1"))
q.fu <- quantile(prodlim(surv_formula, data=clin_data, reverse=T, conf.int = 0.95))

calculate_median_fu(clin_data, "Product")
calculate_median_fu(clin_data, "therapy_before_cart")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Description
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# General description groups
clin_data %>%
  group_by(therapy_before_cart) %>%
  summarize(n = n()) %>%
  mutate(prop=n/sum(n))


# OOS BsAb vs Other BT
clin_data$OOS_binary = factor(clin_data$OOS_binary, levels = c(1, 0))
fisher.test(clin_data$OOS_binary, clin_data$bridging_bin)


# Time from last BsAb infusion to leukapheresis
bite_admin %>%
  filter(days_bite<0) %>%
  group_by(patient_id) %>%
  summarize(n = n()) %>%
  summarize(median = median(n))

final_apheresis = apheresis_oos_df %>%
  group_by(patient_id) %>%
  slice_max(cum_days)

all = apheresis_oos_df %>%
  group_by(patient_id) %>%
  mutate(apheresis = paste0(patient_id, "_apheresis_", row_number()))

bite_admin %>%
  left_join(final_apheresis) %>%
  group_by(patient_id) %>%
  filter(days_bite<cum_days) %>%
  summarize(n = n()) %>%
  summarize(median = median(n))

aphereses_time = bite_admin %>%
  left_join(all) %>%
  group_by(apheresis) %>%
  mutate(diff = days_bite - cum_days) %>%
  filter(days_bite<cum_days) %>%
  slice_min(diff) %>%
  ungroup()

aphereses_time %>%
  wilcox_test(diff ~ oos)

aphereses_time %>%
  summarize(median = median(diff), min=min(diff), max=max(diff))

aphereses_time %>%
  group_by(oos) %>%
  summarize(n=n(), median = median(diff), min=min(diff), max=max(diff))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1B : Swim plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

orr_box <- clin_data %>%
  group_by(`30 day response`) %>%
  summarize(n=n()) %>% 
  mutate(perc=n/sum(n)*100) %>%
  filter(!`30 day response` %in% c("SD/PD")) %>%
  mutate(`30 day response` = factor(`30 day response`, levels=rev(c("CR", "VGPR/PR"))) ) %>%
  ggplot(aes(x="", y=perc, fill=`30 day response`)) +
  geom_bar(position="stack", stat='identity', color="black", alpha=.8) +
  geom_text(aes(label=paste0(`30 day response`, "\n(", round(perc, 0), " %)")), position = position_stack(vjust = 0.5),
            color="black", size=5.5) +
  coord_trans(ylim=c(0,100)) +
  theme_classic() + gtheme(17) + 
  ylab("Day 30 response (%)") + xlab("") +
  theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
  scale_fill_manual("30 day response groups", values=palette_response) +
  guides(fill="none") +
  annotate("text", x=1, y=97, label="ORR: 78%", size=5.5)
orr_box

# Events
death_df <- clin_data[,c("patient", "os", "death")]
colnames(death_df)[2:3] <- c("time", "status") 
death_df$Event <- "Death"

progression_df <- clin_data[,c("patient", "pfs", "status")]
colnames(progression_df)[2] <- c("time") 
progression_df$Event <- "Progression"

event_df <- rbind(death_df, progression_df) %>% filter(status==1) %>%
  group_by(patient) %>% 
  mutate(difference_pr_de = lag(time) - time) %>%
  filter(is.na(difference_pr_de) | difference_pr_de!=0) %>%
  ungroup() %>%
  data.frame()

# Arrows for prolonged remission
arrow_df <- clin_data %>%
  filter(death != 1, progression!= 1) %>%
  dplyr::select(patient, os, Product) %>% 
  mutate(arrow_start = os + 10) %>%
  data.frame()

swim_df <- clin_data[,c("patient", "os", "30 day response", "Product", "refractoriness_group")] %>%
  mutate(response_1 = `30 day response`) %>%
  data.frame()

swim_g <- swimmer_plot(
  df=swim_df,
  id='patient',
  end="os",
  name_fill="Product",
  name_col="Product",
  width=.8,
  base_size=15,
  alpha=.7,
  id_order="response_1") +
  scale_fill_manual(name="CAR-T Product", values=pal_product) +
  swimmer_points(df_points=event_df,id='patient',time='time',name_shape ='Event',size=3.5,fill='white',col='black') + 
  swimmer_arrows(df_arrows=arrow_df,id='patient',arrow_start='arrow_start',type ="open",cex=.65, name_col='Product', 
                 angle = 20, length=0.15, show.legend = F, arrow_positions=c(0.1, 20)) +
  scale_shape_manual(name="Events",values=c(Progression=21, Death=17), breaks=c("Death", "Progression")) +
  scale_color_manual(name="CAR-T Product", values=pal_product) +
  gtheme() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border= element_blank(),
    axis.line.x.bottom = element_line(linewidth=.75),
    axis.title.x = element_text(vjust=-10)
  ) +
  scale_y_continuous(name="Days after CAR T cell infusion", breaks=seq(0, 700, by=100), expand=c(0,1)) +
  annotate("text", x=49.5, y=650, label="Continued response",size=5) +
  annotate("text",x=48, y=650, label=sprintf('\u2192'),size=8)

# Annotation
joined_anno <- clin_data[,c("patient", "therapy_before_cart", "refractoriness_group", "status_before_cart", "30 day response")] %>%
  rename(Bridging = "therapy_before_cart", `Remission before`="status_before_cart", Refractoriness = "refractoriness_group") %>%
  pivot_longer(!patient, names_to="Category", values_to = "Value") %>%
  mutate(
    patient = factor(patient, levels =  levels(swim_g$data$patient)),
    Category = factor(Category, levels=c("Bridging", "Refractoriness", "Remission before", "30 day response"))) %>%
  ggplot(aes(x=patient, y=Category, fill=Value)) +
  geom_tile(height=0.9, width=.85, color="grey90") + 
  coord_flip() +
  theme_minimal(15) +
  gtheme() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=50, hjust=1),
    axis.text.y = element_text(size=13),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.y = element_line(linewidth=.75),
    axis.ticks.y = element_line()
  ) +
  scale_fill_manual(name="Bridging therapy", values=c(anno_nejm, palette_response, refractoriness_col), breaks=names(anno_nejm), labels = c("Bispecific Ab", "anti-CD38", "Chemotherapy", "anti-SLAMF7"), guide = guide_legend(order = 2)) +
  guides(fill=guide_legend(override.aes = list(size=0.5)))  +
  new_scale_fill() +
  new_scale_color() +
  geom_point(aes( shape=NA, color=Value)) +
  scale_color_manual(name="Response", values=c(anno_nejm, palette_response), breaks=names(palette_response)) +
  guides(color=guide_legend(override.aes = list(shape=16,size=9.5)), order = 1) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(aes(color = Value, fill = Value), alpha=0) +
  scale_color_manual(name = "Refractoriness", values=refractoriness_col, breaks=names(refractoriness_col), guide_legend(order = 1)) +
  guides(color = guide_legend(override.aes=list(alpha=1, size=6, shape = 17)), fill = "none", shape="none")

swim_plot_cohort <- (joined_anno | (swim_g + theme(axis.title.x = element_text(vjust=15))) ) +
  plot_layout(widths=c(1.8, 9.5), guides="collect")

swim_plot_cohort_box <- cowplot::ggdraw() +
  cowplot::draw_plot(swim_plot_cohort) +
  cowplot::draw_plot(orr_box, x = 0.55, y = .23, width = .21, height = .36)

ggsave2(
  "figures/main/figure_1/swim_plot.png",
  plot = swim_plot_cohort_box,
  width = 170, height = 162.5, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 2C : Numbers of prior lines of therapy
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tapply(clin_data$n_lines_treatment, clin_data$therapy_before_cart, summary)

lines_treatment_violins = clin_data %>%
  mutate(therapy_before_cart = factor(therapy_before_cart, levels = levels(therapy_before_cart))) %>%
  ggplot(aes(x=therapy_before_cart, y=n_lines_treatment, fill=therapy_before_cart)) +
  geom_violin(alpha = .75, draw_quantiles = .5) +
  geom_jitter(size=2, width=.25, height = 0, alpha=.5) +
  scale_fill_manual(values=anno_nejm) +
  scale_color_manual(values=anno_nejm) +
  guides(fill = "none", color = "none") +
  mytheme_grid(14) +
  scale_y_continuous(breaks = seq(0, 16, 2)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_compare_means(label = "p.format", label.x = .8, label.y = 14, size=4.5) +
  ylab("Number of prior treatment lines") +
  facet_grid(.~paste0("Prior lines of treatment"))

ggsave2(
  "figures/main/figure_1/lines_treatment_violins.png",
  plot = lines_treatment_violins,
  width = 70, height = 100, dpi = 400, bg = "white", units = "mm", scale = 1.6
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1C : Timeline bridging
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

apheresis_oos_df <- bridging.obj$bite_leukapheresis %>%
    mutate(
      patient = as.character(patient),
      oos = factor(oos, levels = c("Successful", "OOS"))
    )

bite_admin <- bridging.obj$bite_admin
bite_treat <- bite_admin %>%
  group_by(patient) %>% 
  filter(days_bite == min(days_bite)) %>% 
  arrange(therapy, days_bite) %>% 
  rename(days_before = "days_bite")

bite_admin$patient <- as.character(bite_admin$patient)
bite_treat$patient <- as.character(bite_treat$patient)

response_groups_str <- "'nCR'='CR'; 'VGPR'='VGPR/PR'; 'PR'='VGPR/PR'; 'SD'='SD/PD'; 'PD'='SD/PD'"

deaths <- clin_data %>%
  filter(patient_id %in% unique(bite_admin$patient_id), death==1) %>% 
  mutate(patient= gsub("Patient0", "", patient_id), event="Death") %>%
  dplyr::select(patient, event, os) %>%
  rename(time="os")

progression <- clin_data %>%
  filter(patient_id %in% unique(bite_admin$patient_id), status==1) %>% 
  mutate(patient= gsub("Patient0", "", patient_id), event="Progression") %>%
  dplyr::select(patient, event, pfs) %>% 
  rename(time="pfs")

asct = data.frame(patient=22, event="ASCT", time=150)

event_df <- rbind(deaths, progression, asct)

before_cart <- clin_data[,c("patient_id", "status_before_cart", "os", "death", "pfs", "status")] %>% 
  filter(patient_id %in% unique(bite_admin$patient_id)) %>%
  mutate(
    patient = gsub("Patient0", "", patient_id),
    start=-5,
    end = case_when(
      os < 30 ~ os,
      pfs < 30 & status == 1 ~ pfs,
      .default = 30
    )
  )
before_cart[which(before_cart$patient_id=="Patient022"),]$end = 42
before_cart[which(before_cart$patient_id=="Patient031"),]$end = 62


response_30 <- clin_data[,c("patient_id", "os", "status", "pfs", "death","30 day response")] %>% 
  filter(patient_id %in% unique(bite_admin$patient_id), os>=30) %>% 
  mutate(patient= gsub("Patient0", "", patient_id)) %>%
  mutate(
    start=ifelse(pfs<30, pfs, 30),
    end=case_when(
      os < 90 ~ os,
      pfs < 90 & pfs > 30 & status==1 ~ pfs,
      .default=90
    )
  )
response_30[which(response_30$patient_id=="Patient022"),]$start = 42
response_30[which(response_30$patient_id=="Patient031"),]$start = 62

response_60 <- clin_data[,c("patient_id", "os", "status", "pfs", "death","3 month response")] %>% 
  filter(patient_id %in% unique(bite_admin$patient_id), os>=90, !is.na(`3 month response`)) %>% 
  mutate(patient= gsub("Patient0", "", patient_id)) %>%
  mutate(
    start=ifelse(pfs<90 & pfs>30 & status==1, pfs, 90),
    end=ifelse(pfs > 90 & status == 1, pfs, os)
  ) %>% 
  mutate(end = ifelse(os>220, 220, os))

response_salvage_one <- clin_data[,c("patient_id", "best_response_first_salvage", "os", "death", "pfs")] %>%
  mutate(patient= gsub("Patient0", "", patient_id)) %>%
  filter(patient_id %in% unique(bite_admin$patient_id), !is.na(best_response_first_salvage), patient!=31) %>% 
  mutate(
    best_response_first_salvage = car::recode(best_response_first_salvage, response_groups_str),
    start = pfs,
    end=150
  ) %>% 
  rename(response = "best_response_first_salvage")

response_salvage_two <- data.frame(
  patient="22", response="VGPR/PR", start=150, end=220
)

arrow_df <- clin_data %>% 
  filter(patient_id %in% unique(bite_admin$patient_id), death != 1) %>%
  mutate(patient= as.character(gsub("Patient0", "", patient_id)), start=ifelse(os>220, 222.5, os+2.5)) %>% 
  dplyr::select(patient, start)

ylabels <- function(breaks){
  labels <- paste0("P", breaks)
  return(labels)
}

unique(bite_treat$patient)

reverse_levels = function(x) {
  x[["patient"]] = factor(x[["patient"]], levels = rev(c("22", "24", "31", "32", "33", "43", "48", "50", "52", "55")))
  x
}

event_df = reverse_levels(event_df)
bite_admin = reverse_levels(bite_admin)
bite_treat = reverse_levels(bite_treat)
before_cart = reverse_levels(before_cart)
response_30 = reverse_levels(response_30)
response_60 = reverse_levels(response_60)
response_salvage_one = reverse_levels(response_salvage_one)

bite_treat$therapy = factor(bite_treat$therapy, levels = c("teclistamab", "talquetamab"))
event_df$event = factor(event_df$event, levels = c("ASCT", "Progression", "Death"))

bite_timeline <- bite_admin %>%
  ggplot(aes(y=patient, yend=patient)) +
  geom_segment(data=before_cart, aes(y=patient, x=start, xend=end, color=status_before_cart), linewidth=8, alpha=.6) +
  geom_segment(data=response_30, aes(y=patient, x=start, xend=end, color=`30 day response`), linewidth=8, alpha=.6) +
  geom_segment(data=response_60, aes(y=patient, x=start, xend=end, color=`3 month response`), linewidth=8, alpha=.6) +
  geom_segment(data=response_salvage_one, aes(y=patient, x=start, xend=end, color=response), linewidth=8, alpha=.6) +
  geom_segment(data=response_salvage_two, aes(y=patient, x=start, xend=end, color=response), linewidth=8, alpha=.6) +
  geom_segment(data=arrow_df, aes(x=start, y=patient,yend=patient, xend=start+10), linewidth=0.9, lineend="round", linejoin="round", arrow=arrow(length=unit(0.125, "inches")), alpha=.75, inherit.aes=F) +
  scale_color_manual(name="Response", values=palette_response) +
  guides(color=guide_legend(nrow=3, order=2, title.position="top", override.aes=list(alpha=.4, shape=17))) +
  ggnewscale::new_scale_color() +
  geom_segment(data=bite_treat, aes(x=days_before, xend=-5), color="grey87", linewidth=8) +
  geom_segment(aes(x=days_bite, xend=days_bite+1.55, color=therapy), linewidth=8) +
  scale_color_manual("BsAb", values=c("#EE7733", "#994455")) +
  geom_vline(xintercept=0, linetype="dashed") +
  theme_classic() +
  gtheme(18) +
  ylab("Patient") +
  scale_x_continuous(lim=c(-270, 240), breaks=c(seq(-300, 250, by=60))) +
  xlab("Days after CAR T cell infusion") +
  scale_y_discrete(expand=c(0,1.5), labels = ylabels) +
  theme(legend.position="bottom", axis.title.y = element_blank(), plot.title=element_text(face="bold")) +
  annotate("text", x=-80, y=11.25, label="Bridging", size=7) +
  annotate("text", x=70, y=11.25, label="Salvage", size=7) +
  annotate("label", x=0, y=11.25, label="CAR-T", size=7, label.size = NA) +
  guides(color=guide_legend(nrow=2, order=1, title.position="top")) +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale_color() +
  geom_star(data=apheresis_oos_df, aes(x=cum_days,y=patient, fill=oos), inherit.aes = F, size=3.5, starshape=1) +
  scale_fill_manual(name = "Leukapheresis",values=c("black", "firebrick1")) +
  guides(fill=guide_legend( title.position="top", nrow=2, order = 3)) +
  geom_point(data=event_df, aes(y=patient, x=time, shape=event), inherit.aes=F, size=3.5, fill="white") +
  scale_shape_manual(name="Event", values=c(Death=17, Progression=21, ASCT=23)) +
  guides( shape=guide_legend(nrow=3, title.position="top"), fill=guide_legend(order=3, nrow=2, title.position="top")) +
  ggtitle("BsAb infusion")

ggsave2(
  "figures/main/figure_1/bite_timeline.png",
  bite_timeline,
  width = 170, height = 130, dpi = 400, bg = "white", units = "mm", scale = 1.5
)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fig. 1D : Median time from first leukapheresis to> infusion
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

tapply(clin_data$time_apheresis_cart, clin_data$therapy_before_cart, summary)

violins_time_apheresis_cart <- clin_data %>%
  mutate(therapy_before_cart = factor(therapy_before_cart, levels = levels(therapy_before_cart))) %>%
  ggplot(aes(x=therapy_before_cart, y=time_apheresis_cart, fill=therapy_before_cart)) +
  geom_violin(alpha = .8, draw_quantiles = .5) +
  geom_jitter(size=1.5, width=.25, alpha=.5) +
  scale_fill_manual(values=anno_nejm) +
  scale_color_manual(values=anno_nejm) +
  guides(fill = "none", color = "none") +
  mytheme_grid(14) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_compare_means(label = "p.format", label.x = .8, label.y = 300, size=4.5) +
  geom_pwc(label = "p.signif", hide.ns = "p") +
  ylab("Days") +
  facet_grid(.~paste0("Time from apheresis to infusion"))

ggsave2(
  "figures/main/figure_1/time_apheresis_cart_violins.png",
  plot = violins_time_apheresis_cart,
  width = 70, height = 100, dpi = 400, bg = "white", units = "mm", scale = 1.6
)