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

# Additional experiments ####
## Comparison of different integration methods (BASS) ####
library(BASS)
library(Matrix)

load("outputs/Domain/retina_st_domain.RData")

C <- 10
R <- 9

retina.list <- SplitObject(retina.st.domain, split.by = "sample.name")
cntm <- list()
xym <- list()
for (retina.i in retina.list) {
  cntm <- c(
    cntm, list(retina.i@assays$Spatial@data[gene.list, ])
  )
  xym <- c(
    xym, list(retina.i@images[[paste0("slice_", retina.i$slice.name[1])]]@coordinates[4:5])
  )
}

set.seed(0)
BASS <- createBASSObject(
  cntm, xym, C = C, R = R,
  beta_method = "SW", init_method = "mclust"
)
BASS <- BASS.preprocess(
  BASS, doLogNormalize = FALSE, doBatchCorrect = TRUE,
  geneSelect = "hvgs", nHVG = length(gene.list), doPCA = TRUE,
  scaleFeature = TRUE, nPC = 30
)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS)

domain <- data.frame(row.names = Cells(retina.st.domain))
domain$domain <- NULL
spot.names <- NULL
for (i in seq_along(BASS@X)) {
  domain[colnames(BASS@X[[i]]), "domain"] <- as.integer(BASS@results$z[[i]])
  spot.names <- c(spot.names, colnames(BASS@X[[i]]))
}
embedding <- t(BASS@X_run)
rownames(embedding) <- spot.names
retina.st.domain$domain <- domain[Cells(retina.st.domain), "domain"]
Idents(retina.st.domain) <- retina.st.domain$domain

retina.st.domain@reductions$bass <- CreateDimReducObject(
  embedding, assay = "Spatial", key = "bass_"
)
retina.st.domain <- RunUMAP(
  retina.st.domain, dims = 1:20, reduction = "bass", reduction.name = "umap.bass"
)

save(
  retina.st.domain,
  file = "outputs/Domain/retina_st_domain_bass.RData"
)

## Comparison of different integration methods (GraphST) ####
load("outputs/Domain/retina_st_domain.RData")

mat <- as.data.frame(retina.st.domain@assays$Spatial@counts[gene.list, ])
write.csv(mat, file = "outputs/Domain/GraphST/data.csv")
batch <- data.frame(
  spot = colnames(retina.st.domain),
  batch = retina.st.domain@meta.data[, "sample.name"]
)
write.csv(batch, file = "outputs/Domain/GraphST/meta.csv")
coo <- list()
for (slice.i in slice.list) {
  obj.i <- subset(retina.st.domain, subset = slice.name == slice.i)
  coo.i <- obj.i@images[[paste0("slice_", slice.i)]]@coordinates
  coo <- c(coo, list(coo.i))
}
coo <- do.call("rbind", coo)[Cells(retina.st.domain), 4:5]
write.csv(coo, file = "outputs/Domain/GraphST/coo.csv")

emb <- read.csv("outputs/Domain/GraphST/emb.csv", sep = " ", header = FALSE)
rownames(emb) <- colnames(retina.st.domain)
retina.st.domain@reductions$graphst <- CreateDimReducObject(
  as.matrix(emb), assay = "Spatial", key = "graphst_"
)
retina.st.domain <- RunUMAP(
  retina.st.domain, dims = 1:20, reduction = "graphst", reduction.name = "umap.graphst"
)

domain <- as.integer(read.csv("outputs/Domain/GraphST/domain.csv", header = FALSE)$V1)
retina.st.domain$domain <- domain
Idents(retina.st.domain) <- retina.st.domain$domain

save(
  retina.st.domain,
  file = "outputs/Domain/retina_st_domain_graphst.RData"
)

## Comparison of different integration methods (Seurat-CCA) ####
load("outputs/Domain/retina_st_domain.RData")
retina.list <- SplitObject(retina.st.domain, split.by = "sample.name")

retina.st.anchors <- FindIntegrationAnchors(
  object.list = retina.list, reduction = "cca",
  anchor.features = gene.list
)
retina.st.inte <- IntegrateData(anchorset = retina.anchors) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
retina.st.inte@reductions$pca@assay.used <- "Spatial"
retina.st.domain@reductions$cca <- retina.st.inte@reductions$pca
retina.st.domain <- RunUMAP(
  retina.st.domain, dims = 1:30, reduction = "cca", reduction.name = "umap.cca"
)
retina.st.domain <- FindNeighbors(
  retina.st.domain, reduction = "cca", dims = 1:30,
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
retina.st.domain@graphs$Spatial_mix_cca <- graph

res.seq <- seq(0.5, 1.5, by = 0.1)

retina.st.domain <- FindClusters(
  retina.st.domain, graph.name = "Spatial_mix_cca",
  resolution = res.seq
)

retina.st.domain$domain <- retina.st.domain$Spatial_mix_cca_res.1.1
Idents(retina.st.domain) <- retina.st.domain$domain

save(
  retina.st.domain,
  file = "outputs/Domain/retina_st_domain_cca.RData"
)

## Comparison of different integration methods (scVI) ####
load("outputs/Domain/retina_st_domain.RData")

mat <- as.data.frame(retina.st.domain@assays$Spatial@counts[gene.list, ])
write.csv(mat, file = "outputs/Domain/scVI/data.csv")
batch <- data.frame(
  spot = colnames(retina.st.domain),
  batch = retina.st.domain@meta.data[, "sample.name"]
)
write.csv(batch, file = "outputs/Domain/scVI/meta.csv")

emb <- read.csv("outputs/Domain/scVI/emb.csv", sep = " ", header = FALSE)
rownames(emb) <- colnames(retina.st.domain)
retina.st.domain@reductions$scvi <- CreateDimReducObject(
  as.matrix(emb), assay = "Spatial", key = "scvi_"
)
retina.st.domain <- RunUMAP(
  retina.st.domain, dims = 1:30, reduction = "scvi", reduction.name = "umap.scvi"
)
retina.st.domain <- FindNeighbors(
  retina.st.domain, reduction = "scvi", dims = 1:30,
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
retina.st.domain@graphs$Spatial_mix_scvi <- graph

res.seq <- seq(0.5, 1.5, by = 0.1)

retina.st.domain <- FindClusters(
  retina.st.domain, graph.name = "Spatial_mix_scvi",
  resolution = res.seq
)

retina.st.domain$domain <- retina.st.domain$Spatial_mix_scvi_res.1.1
Idents(retina.st.domain) <- retina.st.domain$domain

save(
  retina.st.domain,
  file = "outputs/Domain/retina_st_domain_scvi.RData"
)

## iLISI for Harmony / Seurat-CCA / scVI ####
library(lisi)

load("outputs/Domain/retina_st_domain.RData")
meta <- retina.st.domain@meta.data
lisi.harmony <- data.frame(
  iLISI = compute_lisi(
    retina.st.domain@reductions$harmony@cell.embeddings, meta, "sample.name"
  )$sample.name,
  method = "Harmony"
)

load("outputs/Domain/retina_st_domain_bass.RData")
meta <- retina.st.domain@meta.data
lisi.bass <- data.frame(
  iLISI = compute_lisi(
    retina.st.domain@reductions$bass@cell.embeddings, meta, "sample.name"
  )$sample.name,
  method = "BASS"
)

load("outputs/Domain/retina_st_domain_graphst.RData")
meta <- retina.st.domain@meta.data
lisi.graphst <- data.frame(
  iLISI = compute_lisi(
    retina.st.domain@reductions$graphst@cell.embeddings, meta, "sample.name"
  )$sample.name,
  method = "GraphST"
)

load("outputs/Domain/retina_st_domain_cca.RData")
meta <- retina.st.domain@meta.data
lisi.cca <- data.frame(
  iLISI = compute_lisi(
    retina.st.domain@reductions$cca@cell.embeddings, meta, "sample.name"
  )$sample.name,
  method = "Seurat-CCA"
)

load("outputs/Domain/retina_st_domain_scvi.RData")
meta <- retina.st.domain@meta.data
lisi.scvi <- data.frame(
  iLISI = compute_lisi(
    retina.st.domain@reductions$scvi@cell.embeddings, meta, "sample.name"
  )$sample.name,
  method = "scVI"
)

lisi.raw <- data.frame(
  iLISI = compute_lisi(
    t(retina.st.domain@assays$Spatial@scale.data), meta, "sample.name"
  )$sample.name,
  method = "Raw"
)

lisi <- do.call("rbind", list(lisi.harmony, lisi.bass, lisi.graphst, lisi.cca, lisi.scvi, lisi.raw))
save(lisi, file = "outputs/Domain/lisi_6_methods.RData")

## 10-fold cross-validation for Harmony / Seurat-CCA / scVI ####
load("outputs/Domain/retina_st_domain.RData")
expr <- retina.st.domain@assays$Spatial@scale.data
domain.harmony <- retina.st.domain$domain

load("outputs/Domain/retina_st_domain_bass.RData")
domain.bass <- retina.st.domain$domain

load("outputs/Domain/retina_st_domain_graphst.RData")
domain.graphst <- retina.st.domain$domain

load("outputs/Domain/retina_st_domain_cca.RData")
domain.cca <- retina.st.domain$domain

load("outputs/Domain/retina_st_domain_scvi.RData")
domain.scvi <- retina.st.domain$domain

domain <- cbind(domain.harmony, domain.bass, domain.graphst, domain.cca, domain.scvi)

write.csv(expr, file = "outputs/Domain/CrossVal/expr.csv")
write.csv(domain, file = "outputs/Domain/CrossVal/domain.csv")
