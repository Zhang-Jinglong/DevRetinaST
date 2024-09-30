# Packages & settings ####
library(ggplot2)
library(ggpubr)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/ImageProcessing/thickness_anno.RData")

thickness.anno$sample.name <- factor(thickness.anno$sample.name, levels = sample.list)
thickness.anno$layer <- factor(thickness.anno$layer, levels = c("NBL", "GCL"))
thickness.anno$proportion <- thickness.anno$thickness.layer / thickness.anno$thickness.retina

# Retina thickness (Fig. S1a) ####
ggplot(thickness.anno) +
  geom_boxplot(
    aes(x = sample.name, y = thickness.retina, color = sample.name),
    outlier.size = 0.5
  ) +
  geom_hline(
    yintercept = ceiling(median(thickness.anno$thickness.retina)),
    linetype = "dashed", color = "darkred"
  ) +
  scale_y_continuous(
    breaks = c(100, ceiling(median(thickness.anno$thickness.retina)), 200, 300)
  ) +
  scale_color_manual(values = sample.colors) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  "outputs/Visualization/fig-S1a.pdf",
  width = 5.5, height = 4
)

# Layer thickness proportion (Fig. S1b) ####
ggplot(
  thickness.anno,
  aes(x = sample.name, y = proportion, color = layer)
) +
  geom_boxplot(outlier.size = 0.5) +
  stat_compare_means(
    aes(method = "wilcox"),
    label = "p.signif"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  "outputs/Visualization/fig-S1b.pdf",
  width = 5.5, height = 4
)
