# Packages & settings ####
library(Seurat)
library(dplyr)
library(harmony)
library(clustree)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/SVG/svg_list.RData")

# Pre-processing ####
DefaultAssay(retina.st.filter) <- "Spatial"

gene.list <- unique(do.call("rbind", svg.list)$gene.symbol)

retina.st.domain <- NormalizeData(
  retina.st.filter, normalization.method = "LogNormalize",
  scale.factor = 10000, verbose = FALSE
) %>% ScaleData(
  features = gene.list, verbose = FALSE
) %>% RunPCA(
  features = gene.list, npcs = 50, verbose = FALSE
)

# Integration ####
retina.st.domain <- RunHarmony(
  retina.st.domain, group.by.vars = "sample.name",
  reduction.use = "pca", reduction.save = "harmony"
)
retina.st.domain <- RunUMAP(
  retina.st.domain, dims = 1:50, reduction = "harmony"
)

# Build graph ####
retina.st.domain <- FindNeighbors(
  retina.st.domain, reduction = "harmony", dims = 1:50,
  k.param = 15, verbose = FALSE
)

graph <- retina.st.domain@graphs$Spatial_snn
for (slice.i in slice.list) {
  obj.i <- subset(retina.st.domain, subset = slice.name == slice.i)
  coo.i <- obj.i@images[[paste0("slice_", slice.i)]]@coordinates
  coo.i <- coo.i[Cells(obj.i), ]
  
  dist.i <- create.spatial.graph(coo.i)
  graph.sub <- graph[rownames(dist.i), colnames(dist.i)] + dist.i
  graph[rownames(dist.i), colnames(dist.i)] <- graph.sub
}
graph <- graph / 2

graph <- as.Graph(graph)
graph@assay.used <- "Spatial"
retina.st.domain@graphs$Spatial_mix <- graph

# Domain detection ####
res.seq <- seq(0.5, 1.5, by = 0.1)

retina.st.domain <- FindClusters(
  retina.st.domain, graph.name = "Spatial_mix",
  resolution = res.seq
)

# Best resolution
clustree(retina.st.domain, prefix = "Spatial_mix_res.")
diff.res.num.list <- diff.res.num(retina.st.domain, 0.5, 1.5, 0.1)
save(
  diff.res.num.list,
  file = "outputs/Domain/diff_res_num_list.RData"
)

retina.st.domain$domain <- retina.st.domain$Spatial_mix_res.1.1
Idents(retina.st.domain) <- retina.st.domain$domain

# Save data ####
save(
  retina.st.domain,
  file = "outputs/Domain/retina_st_domain.RData"
)
