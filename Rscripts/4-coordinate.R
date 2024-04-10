# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")

# Create new coordinate ####
spatial.coord.new <- data.frame(
  layer = radial.anno[Cells(retina.st.filter), "layer"],
  quantile = radial.anno[Cells(retina.st.filter), "quantile"],
  degree = degree.df[Cells(retina.st.filter), "degree"],
  row.names = Cells(retina.st.filter)
)
nbl.idx <- spatial.coord.new$layer == "NBL"
spatial.coord.new$quantile[nbl.idx] <- 7 - spatial.coord.new$quantile[nbl.idx]
gcl.idx <- spatial.coord.new$layer == "GCL"
spatial.coord.new$quantile[gcl.idx] <- 6 - spatial.coord.new$quantile[gcl.idx]

spatial.coord.new$coord.x <- (
  - spatial.coord.new$quantile * cos(spatial.coord.new$degree / 180 * pi)
)
spatial.coord.new$coord.y <- (
  spatial.coord.new$quantile * sin(spatial.coord.new$degree / 180 * pi)
)
spatial.coord.new <- spatial.coord.new[, 4:5]

# Save data ####
save(
  spatial.coord.new,
  file = "outputs/ImageProcessing/coordinate.RData"
)
