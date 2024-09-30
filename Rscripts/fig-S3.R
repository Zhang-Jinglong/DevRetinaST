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
load("outputs/DEG/sublayer_deg.RData")
load("outputs/DEG/sublayer_stage_deg.RData")
load("outputs/DEG/sublayer_deg_go.RData")

retina.st.domain$radial <- radial.anno[Cells(retina.st.domain), "layer"]

# GCL & NBL marker (Fig. S3a) ####
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
    ) +
    facet_grid(. ~ Gene)
  pic
}

df.1 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["TUBB3",],
  Area = retina.st.domain$radial,
  Gene = "TUBB3"
) # 5.396738e-158
pic.a1 <- marker.plot(df.1, "less")

df.2 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SNCG",],
  Area = retina.st.domain$radial,
  Gene = "SNCG"
) # 3.067698e-167
pic.a2 <- marker.plot(df.2, "less")

df.3 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SOX2",],
  Area = retina.st.domain$radial,
  Gene = "SOX2"
) # 3.024096e-26
pic.a3 <- marker.plot(df.3, "greater")

df.4 <- data.frame(
  Exp = retina.st.domain@assays$Spatial@data["SOX9",],
  Area = retina.st.domain$radial,
  Gene = "SOX9"
) # 9.519671e-30
pic.a4 <- marker.plot(df.4, "greater")

pic.a1 + pic.a2 + pic.a3 + pic.a4 + plot_layout(widths = c(1, 1, 1, 1))

ggsave(
  "outputs/Visualization/fig-S3a.pdf",
  width = 5, height = 3.5
)

# Sublayer DEG (Fig. S3b) ####
deg.num.df <- as.data.frame(table(sublayer.stage.deg[, c("layer", "trans")]))
deg.num.df$layer <- factor(deg.num.df$layer, levels = radial.list)
deg.num.df$trans <- factor(deg.num.df$trans, levels = paste(sample.list[1:5], sample.list[2:6], sep = "."))

ggplot(deg.num.df) +
  geom_bar(
    aes(x = trans, y = Freq + 1, fill = trans),
    stat = "identity", position = "dodge"
  ) +
  facet_grid(~layer) +
  geom_text(
    aes(x = trans, y = Freq + 1, label = Freq, group = trans),
    position = position_dodge(width = 0.9),
    hjust = -0.3, angle = 90, size = 3
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  ylab("Number of DEGs") +
  ylim(0, 200) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(size = 0.1),
    panel.border = element_rect(linewidth = 0.1)
  )

ggsave(
  "outputs/Visualization/fig-S3b.pdf",
  width = 6, height = 3
)

# Sublayer DEG GO (Fig. S3c) ####
dotplot(deg.go, label_format = 50, font.size = 9, showCategory = 5) +
  ggtitle("Gene ontology (Biological process)")

ggsave(
  "outputs/Visualization/fig-S3c.pdf",
  width = 12, height = 6
)
