# Packages & settings ####
rm(list = ls())
source("scripts/utils.R")

# Load GTF file ####
ref.gtf <- read.csv(ref.gtf.file, sep = "\t", row.names = 5)

# Additional experiments ####
## GSE229682 ####
data <- read.csv("data/Bulk/GSE229682/GSE229682_Gene_CPM_MSTR.tsv", sep = "\t", row.names = 1)
data <- data[, 10:24]

data <- data[intersect(rownames(ref.gtf), rownames(data)),]
gene.symbol <- ref.gtf[rownames(data), "gene.symbol"]
rownames(data) <- gene.symbol
retina.bulk <- data

save(retina.bulk, file = "outputs/MergeData/retina_bulk_pub.RData")
