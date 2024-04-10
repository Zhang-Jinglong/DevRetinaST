# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/DEG/nbl_gcl_deg.RData")
load("outputs/DEG/nbl_gcl_deg_go.RData")

retina.st.domain$radial <- radial.anno[Cells(retina.st.domain), "layer"]

# GCL & NBL marker & DEG (Fig. S2a-S2b) ####
## GCL & NBL marker ####
marker.plot <- function(df, alternative) {
  df$Area <- factor(df$Area, levels = c("NBL", "GCL"))
  print(wilcox.test(
    x = df$Exp[df.1$Area == "NBL"],
    y = df$Exp[df.1$Area == "GCL"],
    alternative = alternative
  )$p.value)
  
  pic <- ggplot(
    df,
    aes(x = Area, y = Exp, color = Area)
  ) +
    geom_boxplot(
      size = 0.1,
      outlier.alpha = 0
    ) +
    geom_signif(
      comparisons = list(c("NBL", "GCL")),
      map_signif_level = TRUE, color = "black",
      test = wilcox.test, size = 0.1,
      test.args = alternative,
      y_position = 4
    ) +
    theme_bw() +
    ylab("Gene expression") +
    ylim(-0.1, 4.5) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      strip.background = element_rect(linewidth = 0.1),
      panel.border = element_rect(linewidth = 0.1)
    ) + facet_grid(. ~ Gene)
  pic
}

df.1 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["TUBB3", ],
  Area = retina.st.domain$radial,
  Gene = "TUBB3"
) # 5.396738e-158
pic.a1 <- marker.plot(df.1, "less")

df.2 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SNCG", ],
  Area = retina.st.domain$radial,
  Gene = "SNCG"
) # 3.067698e-167
pic.a2 <- marker.plot(df.2, "less")

df.3 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SOX2", ],
  Area = retina.st.domain$radial,
  Gene = "SOX2"
) # 3.024096e-26
pic.a3 <- marker.plot(df.3, "greater")

df.4 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SOX9", ],
  Area = retina.st.domain$radial,
  Gene = "SOX9"
) # 9.519671e-30
pic.a4 <- marker.plot(df.4, "greater")

## GCL & NBL DEG ####
deg.plot <- function(df) {
  df$sample <- factor(df$sample, levels = sample.list)
  df$Freq <- df$Freq + 0.5
  pic <- ggplot(df) +
    geom_bar(
      aes(x = Var1, y = Freq, fill = sample),
      stat = "identity",
      position = "dodge"
    ) +
    geom_text(
      aes(x = Var1, y = Freq, label = Freq - 0.5, group = sample),
      position = position_dodge(width = 0.9),
      hjust = -0.3, angle = 90, size = 3
    ) +
    scale_fill_manual(values = sample.colors) +
    theme_bw() + ylab("Number of DEGs") + ylim(0, 200) +
    theme(
      axis.text = element_text(color = "black"),
      legend.position = "none",
      axis.title.x = element_blank(),
      strip.background = element_rect(size = 0.1),
      panel.border = element_rect(linewidth = 0.1)
    ) + facet_grid(. ~ group)
  pic
}

deg.nbl.num.df$group <- "NBL"
deg.nbl.num.df$Var1 <- factor(deg.nbl.num.df$Var1, levels = radial.list[1:4])
pic.b1 <- deg.plot(deg.nbl.num.df)

deg.gcl.num.df$group <- "GCL"
deg.gcl.num.df$Var1 <- factor(deg.gcl.num.df$Var1, levels = radial.list[5:6])
pic.b2 <- deg.plot(deg.gcl.num.df)

## Subplots concatenation ####
pic.a1 + pic.a2 + pic.a3 + pic.a4 + pic.b1 + pic.b2 +
  plot_layout(widths = c(1, 1, 1, 1, 3, 1.5))

ggsave(
  "outputs/Visualization/fig-S2ab.pdf",
  width = 10, height = 3.5
)

# GCL & NBL DEG GO (Fig. S3c) ####
dotplot(deg.go, label_format = 50, font.size = 9, showCategory = 20) +
  ggtitle("Gene ontology (Biological process)")

ggsave(
  "outputs/Visualization/fig-S2c.pdf",
  width = 15, height = 6
)
