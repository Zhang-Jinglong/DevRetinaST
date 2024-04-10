# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")

retina.st.domain$radial <- radial.anno[Cells(retina.st.domain), "sublayer"]
retina.st.domain$tangential <- degree.df[Cells(retina.st.domain), "degree"]

# Radial annotation (Fig. S1a) ####
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
    paste0("outputs/Visualization/fig-S1a_", slice.i, ".pdf"),
    width = 2, height = 2
  )
}

# Tangential annotation (Fig. S1b) ####
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
    paste0("outputs/Visualization/fig-S1b_", slice.i, ".pdf"),
    width = 2, height = 2
  )
}
