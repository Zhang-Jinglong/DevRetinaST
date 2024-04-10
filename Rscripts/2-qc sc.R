# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# QC for single-cell retina ####
load("outputs/MergeData/retina_sc.RData")

## Percentage of mitochondrial genes ####
retina.sc[["pct.mt"]] <- PercentageFeatureSet(retina.sc, pattern = "^MT-")

VlnPlot(
  retina.sc,
  features = "pct.mt",
  group.by = "sample.name"
) + ylim(-1, 100) + labs(fill = "Sample") +
  theme(axis.title.x = element_blank())

## Feature scatter for nCount_RNA, nFeature_RNA & pct.mt ####
plot1 <- FeatureScatter(
  retina.sc,
  feature1 = "nCount_RNA", feature2 = "pct.mt"
)
plot2 <- FeatureScatter(
  retina.sc,
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
)
plot1 + plot2

## High quality subset ####
retina.sc.qc <- subset(
  retina.sc,
  subset = pct.mt < 10 &
    nCount_RNA > 1000 &
    nFeature_RNA > 300
)

gene.ava <- (
  (rowSums(retina.sc.qc@assays$RNA@counts) >= 10) &
    (!grepl("^MT-", rownames(retina.sc.qc))) &
    (!grepl("^RP[SL]", rownames(retina.sc.qc)))
)
retina.sc.qc <- subset(
  retina.sc.qc,
  feature = rownames(retina.sc.qc)[gene.ava]
)

save(
  retina.sc.qc,
  file = "outputs/QualityControl/retina_sc_filter.RData"
)

# QC for single-cell ciliary body  ####
load("outputs/MergeData/cb_sc.RData")

## Percentage of mitochondrial genes ####
cb.sc[["pct.mt"]] <- PercentageFeatureSet(cb.sc, pattern = "^MT-")

VlnPlot(
  cb.sc,
  features = "pct.mt"
) + ylim(-1, 100) +
  theme(axis.title.x = element_blank())

## Feature scatter for nCount_RNA, nFeature_RNA & pct.mt ####
plot1 <- FeatureScatter(
  cb.sc,
  feature1 = "nCount_RNA", feature2 = "pct.mt"
)
plot2 <- FeatureScatter(
  cb.sc,
  feature1 = "nCount_RNA", feature2 = "nFeature_RNA"
)
plot1 + plot2

## High quality subset ####
cb.sc.qc <- subset(
  cb.sc,
  subset = pct.mt < 40 &
    nCount_RNA > 500 &
    nCount_RNA < 50000 &
    nFeature_RNA > 300
)

gene.ava <- (
  (rowSums(cb.sc.qc@assays$RNA@counts) >= 10) &
    (!grepl("^MT-", rownames(cb.sc.qc))) &
    (!grepl("^RP[SL]", rownames(cb.sc.qc)))
)
cb.sc.qc <- subset(
  cb.sc.qc,
  feature = rownames(cb.sc.qc)[gene.ava]
)

save(
  cb.sc.qc,
  file = "outputs/QualityControl/cb_sc_filter.RData"
)
