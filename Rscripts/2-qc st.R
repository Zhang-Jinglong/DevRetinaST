# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/MergeData/retina_st.RData")

# Percentage of mitochondrial genes ####
retina.st.qc <- merge(
  x = retina.st[[1]],
  y = as.array(retina.st)[-1],
  project = "retina.st"
)
retina.st.qc$sample.name <- factor(
  retina.st.qc$sample.name, levels = sample.list
)
retina.st.qc$pct.mt <- PercentageFeatureSet(
  retina.st.qc, pattern = "^MT-"
)

VlnPlot(
  retina.st.qc,
  features = "pct.mt",
  group.by = "sample.name"
) + ylim(-1, 15) + labs(fill = "Sample") +
  theme(axis.title.x = element_blank())

# Number of reads ####
VlnPlot(
  retina.st.qc,
  features = "nCount_Spatial",
  group.by = "sample.name"
) + theme(axis.title.x = element_blank()) +
  labs(fill = "Sample")

hist(
  retina.st.qc$nCount_Spatial,
  main = "Histogram of nCount_Spatial",
  breaks = 75, xlab = "nCount_Spatial"
)

# High quality subset ####
retina.st.qc <- subset(
  retina.st.qc,
  subset = pct.mt < 10 &
    nCount_Spatial > 1000
)

gene.ava <- (
  (!grepl("^MT-", rownames(retina.st.qc))) &
    (!grepl("^RP[SL]", rownames(retina.st.qc)))
)
retina.st.qc <- subset(
  retina.st.qc,
  feature = rownames(retina.st.qc)[gene.ava]
)

save(
  retina.st.qc,
  file = "outputs/QualityControl/retina_st_qc.RData"
)

# On retina area ####
retina.st.filter <- subset(
  retina.st.qc,
  subset = is.retina == "Retina"
)

gene.ava <- rowSums(retina.st.filter@assays$Spatial@counts) >= 10
retina.st.filter <- subset(
  retina.st.filter,
  feature = rownames(retina.st.filter)[gene.ava]
)

save(
  retina.st.filter,
  file = "outputs/QualityControl/retina_st_filter.RData"
)

# Additional experiments ####
## GSE234035 ####
load("outputs/MergeData/retina_st_pub.RData")

retina.st.pub$pct.mt <- PercentageFeatureSet(
  retina.st.pub, pattern = "^MT-"
)

retina.st.pub <- subset(
  retina.st.pub,
  subset = pct.mt < 10 &
    nCount_Spatial > 100
)

retina.st.pub <- subset(
  retina.st.pub,
  subset = retina == TRUE
)

gene.ava <- (
  (!grepl("^MT-", rownames(retina.st.pub))) &
    (!grepl("^RP[SL]", rownames(retina.st.pub)))
)
retina.st.pub <- subset(
  retina.st.pub,
  feature = rownames(retina.st.pub)[gene.ava]
)

gene.ava <- rowSums(retina.st.pub@assays$Spatial@counts) >= 10
retina.st.pub <- subset(
  retina.st.pub,
  feature = rownames(retina.st.pub)[gene.ava]
)

save(
  retina.st.pub,
  file = "outputs/QualityControl/retina_st_pub.RData"
)
