# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/QualityControl/retina_sc_filter.RData")
load("outputs/QualityControl/retina_st_pub.RData")

## Remove cell types with too few cells (4 cells) ####
table(data.frame(
  retina.sc.qc$cell.type,
  retina.sc.qc$sample.name
))
retina.sc.filter <- subset(
  retina.sc.qc, subset = cell.type != "MG"
)

# Cell2location ####
## Convert input data ####
sc.count <- data.frame(retina.sc.filter@assays$RNA@counts)
write.csv(
  sc.count,
  file = "outputs/Deconvolution/Cell2location/retina_sc.csv"
)
sc.metadata <- data.frame(
  cell_type = retina.sc.filter$cell.type,
  batch = retina.sc.filter$sample.name
)
write.csv(
  sc.metadata,
  file = "outputs/Deconvolution/Cell2location/retina_sc_meta.csv"
)

st.count <- data.frame(retina.st.filter@assays$Spatial@counts)
write.csv(
  st.count,
  file = "outputs/Deconvolution/Cell2location/retina_st.csv"
)
st.metadata <- data.frame(
  batch = retina.st.filter$slice.name
)
write.csv(
  st.metadata,
  file = "outputs/Deconvolution/Cell2location/retina_st_meta.csv"
)

## python 3-deconvolution.py ####
## outputs: "outputs/Deconvolution/Cell2location/cell2location.csv"

# Correct fraction with nuclear density ####
## Load cell2location results ####
pred.cell2loction <- as.matrix(read.csv(
  "outputs/Deconvolution/Cell2location/cell2location.csv",
  row.names = 1, check.names = FALSE
))
pred.cell2loction <- t(t(pred.cell2loction) / colSums(pred.cell2loction))

## Correct fraction ####
load("outputs/ImageProcessing/nuclear_density.RData")
pred.density <- nuclear.density.df[Cells(retina.st.filter), 1]
correct.cell2loction <- t(t(pred.cell2loction) * pred.density)

## Save data ####
save(
  pred.cell2loction, correct.cell2loction,
  file = "outputs/Deconvolution/decon_results.RData"
)

# Additional experiments ####
## Convert public ST data ####
st.count <- data.frame(retina.st.pub@assays$Spatial@counts)
write.csv(
  st.count,
  file = "outputs/Deconvolution/Cell2location/retina_st_pub.csv"
)
st.metadata <- data.frame(
  batch = retina.st.pub$image.prefix
)
write.csv(
  st.metadata,
  file = "outputs/Deconvolution/Cell2location/retina_st_meta_pub.csv"
)

## python 3-deconvolution pub.py ####
## outputs: "outputs/Deconvolution/Cell2location/cell2location_pub.csv"

## Save data ####
pub.cell2loction <- as.matrix(read.csv(
  "outputs/Deconvolution/Cell2location/cell2location_pub.csv",
  row.names = 1, check.names = FALSE
))
pub.cell2loction <- t(t(pub.cell2loction) / colSums(pub.cell2loction))
colnames(pub.cell2loction) <- sub("^.", "", colnames(pub.cell2loction))

save(
  pub.cell2loction,
  file = "outputs/Deconvolution/decon_results_pub.RData"
)
