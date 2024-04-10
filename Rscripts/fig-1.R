# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")

# UMAP for all spots (Fig. 1c) ####
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2],
  radial = radial.anno[Cells(retina.st.domain), "sublayer"],
  tangential = abs(degree.df[Cells(retina.st.domain), "degree"])
)
umap.data$radial <- factor(umap.data$radial, levels = radial.list)

ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = radial,
    size = tangential
  ), alpha = 0.8, shape = 16) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_manual(values = radial.colors) +
  scale_size_continuous(breaks = c(0, 45, 90, 135), range = c(5, 0.5)) +
  coord_fixed(ratio = 2)

ggsave(
  "outputs/Visualization/fig-1c.pdf",
  width = 7, height = 7
)
