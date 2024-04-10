# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Space Ranger metrics (Table S1) ####
## For all spots ####
metric.list <- list()
for (slice.i in slice.list) {
  metric.list[[slice.i]] <- read.csv(
    paste0("data/ST/", slice.i, "/outs/metrics_summary.csv"),
    check.names = FALSE
  )
}
metric.df <- do.call("rbind", metric.list)
colnames(metric.df)[1] <- "Slice"
metric.df$Slice[1:2] <- c("PCW9 & PCW12 early", "PCW12 late")

write.csv(
  metric.df, file = "outputs/Table/table-S1_1.csv",
  row.names = FALSE
)

## For retina ####
load("outputs/QualityControl/retina_st_filter.RData")
metric.df <- data.frame(
  `Sample` = sample.list,
  `Number of Spots On Retina` = 0,
  `Median Genes per Spot On Retina` = 0,
  `Median UMI Counts per Spot On Retina` = 0,
  `Total Genes Detected On Retina` = 0,
  row.names = sample.list,
  check.names = FALSE
)
for (sample.i in sample.list) {
  obj.sub <- subset(retina.st.filter, subset = sample.name == sample.i)
  metric.df[sample.i, 2] <- ncol(obj.sub)
  metric.df[sample.i, 3] <- median(obj.sub$nFeature_Spatial)
  metric.df[sample.i, 4] <- median(obj.sub$nCount_Spatial)
  metric.df[sample.i, 5] <- sum(rowSums(obj.sub@assays$Spatial@counts > 0) > 0)
}
metric.df$Sample[2:3] <- c("PCW12 early", "PCW12 late")

write.csv(
  metric.df, file = "outputs/Table/table-S1_2.csv",
  row.names = FALSE
)

# GCL & NBL DEG (Table S2) ####
load("outputs/DEG/nbl_gcl_deg.RData")
deg.df <- rbind(deg.nbl.df, deg.gcl.df)[, c(8, 1:7)]
colnames(deg.df) <- c(
  "Sample", "p.value", "avg.log2FC", "pct.1", "pct.2",
  "p.adjust", "Sublayer", "Gene"
)
deg.df$Sample <- factor(deg.df$Sample, levels = sample.list)
deg.df$Sublayer <- factor(deg.df$Sublayer, levels = radial.list)
deg.df <- deg.df[order(deg.df$Sublayer, deg.df$Sample, -deg.df$avg.log2FC), ]

write.csv(
  deg.df, file = "outputs/Table/table-S2.csv",
  row.names = FALSE
)

# GCL & NBL DEG GO (Table S3) ####
load("outputs/DEG/nbl_gcl_deg_go.RData")
deg.go.df <- deg.go@result
colnames(deg.go.df) <- c(
  "GO ID", "Description", "GeneRatio", "BgRatio", "p.value", "p.adjust",
  "q.value", "GeneID", "Count"
)

write.csv(
  deg.go.df, file = "outputs/Table/table-S3.csv",
  row.names = FALSE
)

# Validation of cell type marker (Table S4) ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Deconvolution/decon_results.RData")
pre.marker <- read.csv(
  "outputs/Reference/sc_public_markers.csv",
  check.names = FALSE
)
pre.marker$SCC <- 0
pre.marker$p.value <- 0

diff.gene <- setdiff(pre.marker$Gene, rownames(retina.st.domain))
pre.marker$SCC[pre.marker$Gene %in% diff.gene] <- NA
pre.marker$p.value[pre.marker$Gene %in% diff.gene] <- NA

exp.mat <- retina.st.domain@assays$Spatial@data
frac.mat <- pred.cell2loction[, colnames(exp.mat)]
for (marker.i in 1:nrow(pre.marker)) {
  if (!is.na(pre.marker$SCC[marker.i])) {
    exp.i <- exp.mat[pre.marker$Gene[marker.i], ]
    frac.i <- frac.mat[pre.marker$`Cell type`[marker.i], ]
    
    corr.i <- cor.test(
      x = exp.i, y = frac.i, method = "spearman", exact = TRUE
    )
    pre.marker$SCC[marker.i] <- corr.i$estimate
    pre.marker$p.value[marker.i] <- corr.i$p.value
  }
}
pre.marker$p.adjust <- p.adjust(pre.marker$p.value, method = "BH")

pre.marker$Validated <- FALSE
pre.marker$Validated[!is.na(pre.marker$SCC)][
  pre.marker$p.adjust[!is.na(pre.marker$SCC)] < 0.05 &
    pre.marker$SCC[!is.na(pre.marker$SCC)] > 0
] <- TRUE

write.csv(
  pre.marker, file = "outputs/Table/table-S4.csv",
  row.names = FALSE
)

# SVG (Table S5) ####
load("outputs/SVG/svg_list.RData")
for (sample.i in sample.list) {
  df.i <- svg.list[[sample.i]]
  df.i$sample <- sample.i
  svg.list[[sample.i]] <- df.i
}
svg.df <- do.call("rbind", svg.list)[, c(3, 1, 2)]
colnames(svg.df) <- c("Sample", "Gene", "p.adjust")

write.csv(
  svg.df, file = "outputs/Table/table-S5.csv",
  row.names = FALSE
)

# Domain DEG (Table S6) ####
load("outputs/DEG/domain_deg.RData")
colnames(domain.deg) <- c(
  "p.value", "avg.log2FC", "pct.1", "pct.2",
  "p.adjust", "Domain", "Gene"
)
domain.deg$Domain <- factor(domain.deg$Domain, levels = domain.list)
domain.deg <- domain.deg[order(domain.deg$Domain, -domain.deg$avg.log2FC), ]

write.csv(
  domain.deg, file = "outputs/Table/table-S6.csv",
  row.names = FALSE
)

# Domain DEG GO (Table S7) ####
load("outputs/DEG/domain_deg_go.RData")
deg.go.df <- deg.go@compareClusterResult[, -2][, c(2, 1, 3:10)]
colnames(deg.go.df) <- c(
  "GO ID", "Domain", "Description", "GeneRatio", "BgRatio",
  "p.value", "p.adjust", "q.value", "GeneID", "Count"
)

write.csv(
  deg.go.df, file = "outputs/Table/table-S7.csv",
  row.names = FALSE
)

# Gene module (Table S8) ####
me.list <- list()
for (me.i in paste0("ME", 1:4)) {
  gene.i <- as.character(read.csv(
    paste0("outputs/WGCNA/", me.i, ".txt"), header = FALSE
  )$V1)
  
  if (me.i != "ME1") {
    gene.i <- c(gene.i, rep("", length(me.list[["ME1"]]) - length(gene.i)))
  }
  me.list[[me.i]] <- gene.i
}
me.df <- as.data.frame(me.list)
colnames(me.df) <- paste0("Module ", 1:4)

write.csv(
  me.df, file = "outputs/Table/table-S8.csv",
  row.names = FALSE
)

# Gene module GO (Table S9) ####
load("outputs/DEG/me_go.RData")
me.go.df <- me.go@compareClusterResult[, -2][, c(2, 1, 3:10)]
colnames(me.go.df) <- c(
  "GO ID", "Module", "Description", "GeneRatio", "BgRatio",
  "p.value", "p.adjust", "q.value", "GeneID", "Count"
)

write.csv(
  me.go.df, file = "outputs/Table/table-S9.csv",
  row.names = FALSE
)

# Regulatory gene (Table S10) ####
load("outputs/Development/reg_gene_df.RData")
reg.gene.df <- reg.gene.df[, c(8, 9, 7, 1:6, 10)]
colnames(reg.gene.df) <- c(
  "Source sample", "Target sample", "Domain", "p.value", "avg.log2FC",
  "pct.1", "pct.2", "p.adjust", "Gene", "Gene type"
)

reg.gene.df$`Source sample` <- factor(
  reg.gene.df$`Source sample`, levels = sample.list
)
reg.gene.df$Domain <- factor(reg.gene.df$Domain, levels = domain.list)
order.idx <- order(
  reg.gene.df$`Source sample`, reg.gene.df$Domain, -reg.gene.df$avg.log2FC
)
reg.gene.df <- reg.gene.df[order.idx, ]

write.csv(
  reg.gene.df, file = "outputs/Table/table-S10.csv",
  row.names = FALSE
)

# Regulatory gene GO (Table S11) ####
load("outputs/DEG/reg_gene_go.RData")
reg.go.df <- reg.go@result
colnames(reg.go.df) <- c(
  "GO ID", "Description", "GeneRatio", "BgRatio",
  "p.value", "p.adjust", "q.value", "GeneID", "Count"
)

write.csv(
  reg.go.df, file = "outputs/Table/table-S11.csv",
  row.names = FALSE
)

# Disease gene (Table S12) ####
retnet.table <- readxl::read_excel(
  "outputs/Reference/RetNet_20221007.xlsx",
  sheet = "Merge"
)
colnames(retnet.table)[3] <- "Gene"

write.csv(
  retnet.table, file = "outputs/Table/table-S12.csv",
  row.names = FALSE
)

# Spatial communication score (Table S13) ####
lr.table <- read.csv(
  "outputs/CCC/COMMOT/retina_db_filter.csv", row.names = 1
)
rownames(lr.table) <- paste(lr.table$X0, lr.table$X1, sep = "-")

lr.summary.list <- list()
for (sample.i in sample.list) {
  lr.summary.i <- read.csv(
    paste0("outputs/CCC/COMMOT/commot_", sample.i, "_lr.csv"),
    header = FALSE
  )
  colnames(lr.summary.i) <- c("L-R pair", "SCS")
  lr.summary.i$Sample <- sample.i
  
  lr.summary.list[[sample.i]] <- lr.summary.i[, c(3, 1, 2)]
}
lr.summary <- do.call("rbind", lr.summary.list)
lr.summary$`L-R type` <- lr.table[lr.summary$`L-R pair`, "X3"]

lr.summary$Sample <- factor(lr.summary$Sample, levels = sample.list)
lr.summary <- lr.summary[order(lr.summary$Sample, -lr.summary$SCS), ]

write.csv(
  lr.summary, file = "outputs/Table/table-S13.csv",
  row.names = FALSE
)

# CCC (Table S14) ####
cpdb.means <- read.csv(
  "outputs/CCC/CellPhoneDB/cpdb_filter_means.csv",
  row.names = 1
)
rownames(cpdb.means) <- paste(cpdb.means$gene_a, cpdb.means$gene_b, sep = "-")
cpdb.means <- cpdb.means[, c(-1, -2)]

cpdb.pvals <- read.csv(
  "outputs/CCC/CellPhoneDB/cpdb_filter_pvalues.csv",
  row.names = 1
)
rownames(cpdb.pvals) <- paste(cpdb.pvals$gene_a, cpdb.pvals$gene_b, sep = "-")
cpdb.pvals <- cpdb.pvals[, c(-1, -2)][, colnames(cpdb.means)]

df.cols <- colnames(cpdb.means)
df.cols <- gsub("AC.HC.Precur", "AC-HC Precur", df.cols)
df.cols <- gsub("BC.Photo.Precur", "BC-Photo Precur", df.cols)
df.cols <- gsub("[.]", "|", df.cols)
colnames(cpdb.means) <- df.cols
colnames(cpdb.pvals) <- df.cols

cpdb.ava.idx <- colSums(cpdb.pvals < 0.05) > 0
cpdb.means <- cpdb.means[, cpdb.ava.idx]
cpdb.pvals <- cpdb.pvals[, cpdb.ava.idx]

cpdb.cat <- cpdb.means
for (col.i in colnames(cpdb.means)) {
  cpdb.cat[[col.i]] <- paste0(
    cpdb.means[[col.i]], " (",
    cpdb.pvals[[col.i]], ")"
  )
}

write.csv(
  cpdb.cat, file = "outputs/Table/table-S14.csv"
)
