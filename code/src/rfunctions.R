.cran_packages = c("survival", "survminer", "tidyverse", "stats", "ggsci", "prodlim")

## Requiring packages
for (pack in .cran_packages) {
  suppressMessages(require(
    pack,
    quietly = TRUE,
    character.only = TRUE
  ))
}

source("code/src/ggstyles.R")

# pvalue to code
pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*", ""))))
}

calculate_median_fu <- function(df, group) {
  formula <- as.formula(paste0("Surv(os_months, death)~", group))
  q.fu <- quantile(prodlim(formula, data=df, reverse=T, conf.int = 0.95))
  df <- do.call(cbind, lapply(q.fu,  function(x) { df <- data.frame(x) }))
  names(df) <- names(q.fu)
  return(df)
}

km_plot <- function(x, time, status, groups, palette, title, pval_coord, limits=c(0,21), conf = F, break.by = 3, xtitle="Months after CAR-T cell infusion") {
  surv <- as.formula(paste0("Surv(", time, ",", status, ") ~", groups))
  surv_fit <- survfit(surv, data=x)
  surv_fit$call$formula <- surv
  
  print(surv_fit)

  g_surv <- ggsurvplot(surv_fit, data=x, 
                   color= groups, 
                   legend.title=title,
                   size=1,
                   conf.int= conf,
                   xlim=limits,
                   break.time.by=break.by, 
                   pval=T,  pval.coord=pval_coord, pval.size=8,
                   surv.median.line="hv",
                   palette=palette,
                   risk.table=TRUE,
                   tables.col=groups, tables.y.text=FALSE,
                   risk.table.title="No. at risk",risk.table.fontsize=6.5,
                   ggtheme = theme_classic(20),
  )

  g <- g_surv$plot +
    theme_classic() +
    gtheme(18) +
    theme(legend.key.width = unit(1.1,"cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_rect(linewidth=1, fill=NA),
          strip.text.x = element_text(size=18)
        ) + 
      ylab("PFS Probability") + 
      guides(color = guide_legend(override.aes = list(shape = NA))) +
      xlab(xtitle) +
      scale_color_manual(values=palette, limits=levels(x[[groups]]))

  gt <- g_surv$table +
    theme(
        axis.line=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(size = 18, vjust=-2),
        legend.position="none",
    ) +
    xlab("") +
    ylab("") +
    scale_color_manual(values=palette, limits=levels(x[[groups]]))


    g_g <- ggarrange(g + theme(legend.position="top"), gt +theme(plot.margin = unit(c(-0.5, 0, 0, 0), "cm")),
                    ncol=1, nrow=2, heights=c(5,1.4))
}