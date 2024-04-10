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
