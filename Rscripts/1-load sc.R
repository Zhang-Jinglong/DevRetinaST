# Packages & settings ####
library(Seurat)
library(Matrix)

rm(list = ls())
source("scripts/utils.R")

# Load GTF file ####
ref.gtf <- read.csv(ref.gtf.file, sep = "\t", row.names = 5)

# GSE138002 ####
## Load data ####
sc.data.dir <- "data/SC/GSE138002/"
retina.mat <- readMM(paste0(
  sc.data.dir, "GSE138002_Final_matrix.mtx"
))
gene.list <- read.csv(paste0(
  sc.data.dir, "GSE138002_genes.csv"
), sep = ";", row.names = 1)
sc.metadata <- read.csv(paste0(
  sc.data.dir, "GSE138002_Final_barcodes.csv"
), sep = ";", row.names = 1)

## Adjust format ####
rownames(sc.metadata) <- sc.metadata$barcode
colnames(sc.metadata)[c(2, 6, 10)] <- c(
  "sample.name", "sample.type", "cell.type"
)

rownames(retina.mat) <- gene.list$id
retina.mat <- retina.mat[rownames(retina.mat) %in% rownames(ref.gtf), ]
gene.symbol <- ref.gtf[rownames(retina.mat), "gene.symbol"]
rownames(retina.mat) <- gene.symbol

colnames(retina.mat) <- sc.metadata$barcode

## Convert to Seurat object ####
retina.sc <- CreateSeuratObject(
  retina.mat, project = "retina.sc", assay = "RNA",
  meta.data = sc.metadata[, c(2, 6, 10)]
)
map <- setNames(
  c(
    "AC-HC Precur", "AC", "BC-Photo Precur", "BC", "Cone",
    "HC", "MG", 'Neurogenic', "RGC", "Rod", "RPC"
  ),
  c(
    "AC/HC_Precurs", "Amacrine Cells", "BC/Photo_Precurs",
    "Bipolar Cells", "Cones", "Horizontal Cells", "Muller Glia",
    "Neurogenic Cells", "Retinal Ganglion Cells", "Rods", "RPCs"
  )
)
retina.sc$cell.type <- map[retina.sc$cell.type]

## Select retinal samples that meet the time period ####
sc.sample.list <- paste0("Hgw", c(9, 11:17))
retina.sc <- subset(
  retina.sc,
  subset = sample.type == "Whole Retina" & sample.name %in% sc.sample.list
)
retina.sc$sample.name <- paste0(
  "GW", gsub("Hgw", "", retina.sc$sample.name)
)
print("Single-cell data loaded!")
table(data.frame(retina.sc$sample.name, retina.sc$cell.type))

## Save data ####
save(retina.sc, file = "outputs/MergeData/retina_sc.RData")

# GSE206026 ####
## Load data ####
cb.mat <- Read10X(data.dir = "data/SC/GSE206026", gene.column=1)

## Adjust format ####
cb.mat <- cb.mat[rownames(cb.mat) %in% rownames(ref.gtf), ]
gene.symbol <- ref.gtf[rownames(cb.mat), "gene.symbol"]
rownames(cb.mat) <- gene.symbol

## Convert to Seurat object ####
cb.sc <- CreateSeuratObject(counts = cb.mat)

## Save data ####
save(cb.sc, file = "outputs/MergeData/cb_sc.RData")
