# Packages & settings ####
library(Seurat)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# GCL & NBL DEG ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
retina.st.domain$layer <- radial.anno[Cells(retina.st.domain), "layer"]
retina.st.domain$sublayer <- radial.anno[Cells(retina.st.domain), "sublayer"]

DefaultAssay(retina.st.domain) <- "Spatial"
deg.nbl.list <- list()
deg.gcl.list <- list()
deg.nbl.num.list <- list()
deg.gcl.num.list <- list()
for (sample.i in sample.list) {
  ## GCL DEG
  obj.gcl.i <- subset(
    retina.st.domain,
    subset = (sample.name == sample.i) &
      (layer == "GCL")
  )
  Idents(obj.gcl.i) <- obj.gcl.i$sublayer
  deg.gcl.i <- FindAllMarkers(
    obj.gcl.i, test.use = "wilcox", only.pos = TRUE,
    logfc.threshold = 0.25, min.pct = 0.1
  )
  deg.gcl.i <- deg.gcl.i[deg.gcl.i$p_val_adj < 0.05, ]
  
  if (nrow(deg.gcl.i) > 0) {
    deg.gcl.i$sample <- sample.i
    deg.gcl.list[[sample.i]] <- deg.gcl.i
  }
  deg.gcl.i.count <- as.data.frame(table(deg.gcl.i$cluster))
  deg.gcl.i.count$sample <- sample.i
  deg.gcl.num.list[[sample.i]] <- deg.gcl.i.count
  
  ## NBL DEG
  obj.nbl.i <- subset(
    retina.st.domain,
    subset = (sample.name == sample.i) &
      (layer == "NBL")
  )
  Idents(obj.nbl.i) <- obj.nbl.i$sublayer
  deg.nbl.i <- FindAllMarkers(
    obj.nbl.i, test.use = "wilcox", only.pos = TRUE,
    logfc.threshold = 0.25, min.pct = 0.1
  )
  deg.nbl.i <- deg.nbl.i[deg.nbl.i$p_val_adj < 0.05, ]
  
  if (nrow(deg.nbl.i) > 0) {
    deg.nbl.i$sample <- sample.i
    deg.nbl.list[[sample.i]] <- deg.nbl.i
  }
  deg.nbl.i.count <- as.data.frame(table(deg.nbl.i$cluster))
  deg.nbl.i.count$sample <- sample.i
  deg.nbl.num.list[[sample.i]] <- deg.nbl.i.count
}

deg.gcl.df <- do.call("rbind", deg.gcl.list)
deg.nbl.df <- do.call("rbind", deg.nbl.list)
deg.nbl.num.df <- do.call("rbind", deg.nbl.num.list)
deg.gcl.num.df <- do.call("rbind", deg.gcl.num.list)
save(
  deg.gcl.df, deg.nbl.df, deg.gcl.num.df, deg.nbl.num.df,
  file = "outputs/DEG/nbl_gcl_deg.RData"
)

# GCL & NBL DEG GO ####
deg.list <- unique(c(deg.gcl.df$gene, deg.nbl.df$gene))
deg.entrez <- bitr(deg.list, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
deg.go <- enrichGO(
  deg.entrez$ENTREZID, "org.Hs.eg.db", ont = "BP", pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
)

save(
  deg.go,
  file = "outputs/DEG/nbl_gcl_deg_go.RData"
)

# Domain-specific DEG ####
load("outputs/Domain/retina_st_domain.RData")
retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

Idents(retina.st.domain) <- retina.st.domain$domain
domain.deg <- FindAllMarkers(
  retina.st.domain, test.use = "wilcox",
  logfc.threshold = 0.25, min.pct = 0.1
)
domain.deg <- domain.deg[domain.deg$p_val_adj < 0.05, ]

save(
  domain.deg,
  file = "outputs/DEG/domain_deg.RData"
)

# Domain-specific DEG GO ####
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0, ]
deg.list <- unique(domain.deg$gene)
group <- data.frame(
  gene = domain.deg$gene,
  group = domain.deg$cluster
)
deg.entrez <- bitr(deg.list, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
data <- merge(deg.entrez, group, by.x = "SYMBOL", by.y = "gene")

deg.go <- compareCluster(
  ENTREZID ~ group, 
  data = data, 
  fun = "enrichGO", 
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

save(
  deg.go,
  file = "outputs/DEG/domain_deg_go.RData"
)

# Consensus DEG between D4 & D8 ####
load("outputs/DEG/domain_deg.RData")
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0, ]
single.list <- list(
  D4 = unique(domain.deg$gene[domain.deg$cluster == "D4"]),
  D8 = unique(domain.deg$gene[domain.deg$cluster == "D8"])
)

co.deg.d4.d8 <- intersect(single.list$D4, single.list$D8)

save(
  co.deg.d4.d8,
  file = "outputs/DEG/consensus_deg_d4_d8.RData"
)

# Gene module GO ####
me.list <- c()
for (me.i in paste0("ME", 1:4)) {
  gene.i <- read.table(paste0("outputs/WGCNA/", me.i, ".txt"))
  gene.i$group <- me.i
  me.list[[me.i]] <- gene.i
}

me.df <- do.call("rbind", me.list)
me.list <- unique(me.df$V1)
deg.entrez <- bitr(me.list, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
data <- merge(deg.entrez, me.df, by.x = "SYMBOL", by.y = "V1")

me.go <- compareCluster(
  ENTREZID ~ group, 
  data = data, 
  fun = "enrichGO", 
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

save(
  me.go,
  file = "outputs/DEG/me_go.RData"
)

# Reg. gene GO ####
load("outputs/Development/reg_gene_df.RData")

reg.list <- unique(reg.gene.df$gene)
reg.entrez <- bitr(reg.list, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
reg.go <- enrichGO(
  reg.entrez$ENTREZID, "org.Hs.eg.db", ont = "BP", pAdjustMethod = "bonferroni",
  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
)

save(
  reg.go,
  file = "outputs/DEG/reg_gene_go.RData"
)
