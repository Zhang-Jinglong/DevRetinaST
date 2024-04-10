# Packages & settings ####
library(Seurat)
library(dplyr)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_sc_filter.RData")
load("outputs/QualityControl/cb_sc_filter.RData")

# Pre-processing for retina single-cell data ####
retina.sc.filter <- subset(
  retina.sc.qc, subset = cell.type != "MG"
)
DefaultAssay(retina.sc.filter) <- "RNA"

retina.sc.filter <- NormalizeData(
  retina.sc.filter, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = TRUE) %>%
  ScaleData(
    vars.to.regress = c("sample.name", "nCount_RNA"),
    verbose = FALSE
  ) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, seed.use = 0)

save(
  retina.sc.filter,
  file = "outputs/SC/retina_sc_filter.RData"
)

# Pre-processing for ciliary body single-cell data ####
cb.sc.filter <- NormalizeData(
  cb.sc.qc, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
) %>%
  FindVariableFeatures(nfeatures = 2000, verbose = TRUE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, seed.use = 0) %>%
  FindClusters(resolution = 0.5)

marker.list <- c(
  "OPTC", "BEST2", # NPE
  "SLC4A4", "TFPI2", # PE
  "MLANA", "MITF", # Melanocyte
  "MYH11", "ACTA2", # SMC
  "VWF", "PECAM1", "ERG", # EC
  "C1QA", "CD163", # Macrophage
  "PDGFRA", "COL1A1", # Fibroblast
  "LGI4", "CDH19", # Schwann
  "CD3E", "CD2", # T cell
  "MS4A1", "CD79A", # B cell
  "MZB1", # Plasma cell
  "LYZ", # Monocyte
  "S100A8", "S100A9", # Neutrophil
  "MS4A2", "KIT" # Mast cell
)
# VlnPlot(
#   cb.sc.filter, features = marker.list[1:11]
# )
# VlnPlot(
#   cb.sc.filter, features = marker.list[12:19]
# )
# VlnPlot(
#   cb.sc.filter, features = marker.list[20:27]
# )

rename.list <- c(
  `4` = "NPE", `10` = "NPE", `17` = "NPE",
  `16` = "PE", `6` = "Melanocyte", `20` = "Melanocyte",
  `15` = "SMC", `18` = "EC", `1` = "Macrophage",
  `2` = "Macrophage", `3` = "Macrophage", `12` = "Macrophage",
  `13` = "Macrophage",`9` = "Fibroblast", `7` = "Schwann",
  `0` = "T cell", `11` = "T cell", `5` = "B cell",
  `8` = "Monocyte", `14` = "Neutrophil", `19` = "Mast cell"
)
cb.sc.filter$cell.type <- rename.list[
  as.character(cb.sc.filter$seurat_clusters)
]
Idents(cb.sc.filter) <- cb.sc.filter$cell.type

save(
  cb.sc.filter,
  file = "outputs/SC/cb_sc_filter.RData"
)
