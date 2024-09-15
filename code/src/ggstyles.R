.cran_packages = c("tidyverse", "MetBrewer", "RColorBrewer", "ggthemes", "scales",
                   "ggsci", "scico", "ggstar", "patchwork")

## Requiring packages
for (pack in .cran_packages) {
  suppressMessages(require(
    pack,
    quietly = TRUE,
    character.only = TRUE
  ))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# mytheme (credits: Michael Rade @FraunhoferIZI)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

mytheme_grid = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# gtheme
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gtheme <- function(size=18) {
    theme(text=element_text(size=size,colour="black"),
          axis.text=element_text(size=size,colour="black"),
          axis.title=element_text(size=size,colour="black"),
          strip.text.x = element_text(size=size),
          legend.text=element_text(size=size))
}

colors_use.10 = ggthemes::tableau_color_pal("Tableau 10")(10)

# Palette Response groups
palette_response = c(
  "CR" = "#6699CC",
  "VGPR/PR" = "#EECC66",
  "SD/PD" = "#997700"
)

# Palette type of bridging therapy
nejm <- pal_nejm("default")(8)
anno_nejm <- c(
    "Bispecific Ab" = nejm[1],
    "CD38" = nejm[2],
    "Chemotherapy" = nejm[3],
    "SLAMF7" = nejm[4]
)

# refractoriness palette
col.refr = scico(7, palette="vik")[3:7]

# Treatment refractoriness groups
refractoriness_col <- c(
  TCexposed = col.refr[1],
  TCRRMM = col.refr[2],
  PentaRRMM = col.refr[3]
)

pal_product <- c("Ide-cel" = "grey60", "Cilta-cel" = "#004488")

# CRS
pal_crs <- c(
    "no CRS" = demuth[10],
    "CRS" = demuth[2],
    "CRS-toci" = demuth[4],
    "CRS+toci" = demuth[2]
)

# CLonotype groups
clono.col = c(`Hyperexpanded (100 < X <= 1905)` = "#FF4B20", `Large (20 < X <= 100)`="#FFb433",
              `Medium (5 < X <= 20)`="#C6FDEC", `Small (1 < X <= 5)`="#7AC5FF",
              `Single (0 < X <= 1)`="#0348A6")
