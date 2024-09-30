# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load GTF file ####
ref.gtf <- read.csv(ref.gtf.file, sep = "\t", row.names = 5)

st.data.dirs <- paste0("data/ST/", slice.list, "/outs")
retina.st <- list()
for (data.i in 1:length(st.data.dirs)) {
  # Load count data ####
  obj <- suppressWarnings(Load10X_Spatial(
    data.dir = st.data.dirs[data.i],
    use.names = FALSE,
    slice = paste0("slice_", slice.list[data.i])
  ))
  
  # Load metadata ####
  obj$slice.name <- slice.list[data.i]
  
  sample.name <- read.csv(
    paste0(st.data.dirs[data.i], "/annotation/SampleID.csv")
  )
  obj$sample.name <- NA
  obj$sample.name[sample.name$Barcode] <- sample.name$SampleID
  
  is.retina <- read.csv(
    paste0(st.data.dirs[data.i], "/annotation/Retina.csv")
  )
  obj$is.retina <- NA
  obj$is.retina[is.retina$Barcode] <- is.retina$Retina
  
  obj <- RenameCells(
    obj,
    new.names = paste0(slice.list[data.i], ".", Cells(obj))
  )
  
  # Filtering spots & genes ####
  obj <- subset(
    obj,
    cells = Cells(obj)[which(obj$sample.name != "")],
    features = rownames(obj)[rownames(obj) %in% rownames(ref.gtf)]
  )
  
  gene.symbol <- ref.gtf[rownames(obj), "gene.symbol"]
  rownames(obj@assays$Spatial@counts) <- gene.symbol
  rownames(obj@assays$Spatial@data) <- gene.symbol
  rownames(obj@assays$Spatial@meta.features) <- gene.symbol
  
  retina.st[[slice.list[data.i]]] <- obj
  print(paste0(slice.list[data.i], " spatial data loaded!"))
}

# Save data ####
save(retina.st, file = "outputs/MergeData/retina_st.RData")

# Additional experiments ####
## GSE234035 ####
img.dir.list <- list.files("data/ST/GSE234035/GSE234035_RAW")
img.dir.list <- as.data.frame(t(sapply(strsplit(img.dir.list, "_"), "[", 1:4)))
img.dir.list <- unique(img.dir.list)
img.dir.list$prefix <- paste(img.dir.list$V1, img.dir.list$V2, img.dir.list$V3, img.dir.list$V4, sep = "_")
rownames(img.dir.list) <- paste(img.dir.list$V3, img.dir.list$V4, sep = "_")

ind.list <- c("14749_13PCW", "15685_11PCW", "15817_8PCW", "15870_10PCW")
info.list <- c("Barcode", "stage", "section", "seurat_clusters")
count.list <- list()
meta.list <- list()
for (ind.i in ind.list) {
  count.i <- read.csv(paste0(
    "data/ST/GSE234035/GSE234035_", ind.i, "_counts.csv"
  ), row.names = 1)
  colnames(count.i) <- paste(ind.i, colnames(count.i), sep = "_")

  meta.i <- read.csv(paste0(
    "data/ST/GSE234035/GSE234035_", ind.i, "_md.csv"
  ), row.names = 1)
  meta.i <- meta.i[, info.list]
  rownames(meta.i) <- paste(ind.i, rownames(meta.i), sep = "_")
  meta.i$section <- substr(meta.i$section, 1, 1)
  meta.i$image.prefix <- img.dir.list[paste(meta.i$stage, meta.i$section, sep = "_"), "prefix"]

  count.list <- c(count.list, list(count.i))
  meta.list <- c(meta.list, list(meta.i))
}

count <- do.call("cbind", count.list)
meta <- do.call("rbind", meta.list)

count <- count[rownames(count) %in% rownames(ref.gtf), ]
gene.symbol <- ref.gtf[rownames(count), "gene.symbol"]
rownames(count) <- gene.symbol

retina.area <- read.csv("data/ST/GSE234035/retina_area.csv", row.names = 1)
meta$x <- retina.area[rownames(meta), "Image_X"]
meta$y <- retina.area[rownames(meta), "Image_Y"]
meta$retina <- retina.area[rownames(meta), "retina"] == "True"

retina.st.pub <- CreateSeuratObject(counts = count, meta.data = meta, assay = "Spatial")

save(retina.st.pub, file = "outputs/MergeData/retina_st_pub.RData")
