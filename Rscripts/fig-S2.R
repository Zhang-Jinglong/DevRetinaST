# Packages & settings ####
library(Seurat)
library(ggplot2)
library(patchwork)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")

retina.st.domain$radial <- radial.anno[Cells(retina.st.domain), "sublayer"]
retina.st.domain$tangential <- degree.df[Cells(retina.st.domain), "degree"]

# Radial annotation (Fig. S2a) ####
for (slice.i in slice.list) {
  obj <- subset(
    retina.st.domain,
    subset = slice.name == slice.i
  )
  Idents(obj) <- obj$radial
  SpatialDimPlot(
    obj, pt.size.factor = 1.0, stroke = 0.05,
    images = paste0("slice_", slice.i),
    image.alpha = 0.5, crop = FALSE
  ) + labs(fill = "Radial") +
    scale_fill_manual(values = radial.colors) +
    NoLegend()

  ggsave(
    paste0("outputs/Visualization/fig-S2a_", slice.i, ".pdf"),
    width = 2, height = 2
  )
}

# Tangential annotation (Fig. S2b) ####
for (slice.i in slice.list) {
  obj <- subset(
    retina.st.domain,
    subset = slice.name == slice.i
  )
  SpatialFeaturePlot(
    obj, features = "tangential", crop = FALSE,
    images = paste0("slice_", slice.i),
    max.cutoff = 135, min.cutoff = -135,
    pt.size.factor = 1.0, image.alpha = 0.5, stroke = 0.05
  ) + labs(fill = "Tangential") +
    NoLegend()

  ggsave(
    paste0("outputs/Visualization/fig-S2b_", slice.i, ".pdf"),
    width = 2, height = 2
  )
}

# UMAP colored by radial sublayer, tangential angle & timeline (Fig. S2c) ####
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2],
  sample = retina.st.domain$sample.name,
  radial = retina.st.domain$radial,
  tangential = abs(retina.st.domain$tangential)
)
umap.data$radial <- factor(umap.data$radial, levels = radial.list)
umap.data$sample <- factor(umap.data$sample, levels = sample.list)

pic.radial <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = radial
  ), shape = 16) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_manual(values = radial.colors) +
  coord_fixed(ratio = 2)

pic.tangential <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = tangential
  ), shape = 16) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_continuous(breaks = c(0, 45, 90, 135)) +
  coord_fixed(ratio = 2)

pic.sample <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = sample
  ), shape = 16) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_manual(values = sample.colors) +
  coord_fixed(ratio = 2)

(pic.radial | pic.tangential | pic.sample) + plot_layout(widths = c(1, 1, 1))

ggsave(
  "outputs/Visualization/fig-S2c.pdf",
  width = 15, height = 5
)
