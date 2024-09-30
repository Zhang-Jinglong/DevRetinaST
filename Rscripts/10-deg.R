# Packages & settings ####
library(Seurat)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# Sublayer DEG across stages ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
retina.st.domain$sublayer <- radial.anno[Cells(retina.st.domain), "sublayer"]

deg.list <- list()
for (sublayer.i in radial.list) {
  obj.i <- retina.st.domain[, retina.st.domain$sublayer == sublayer.i]

  Idents(obj.i) <- factor(obj.i$sample.name, levels = sample.list)
  deg.list.i <- list()
  for (sample.i in 1:5) {
    deg.i <- FindMarkers(
      obj.i,
      ident.1 = sample.list[sample.i], ident.2 = sample.list[sample.i + 1],
      test.use = "wilcox", only.pos = TRUE,
      logfc.threshold = 0.25, min.pct = 0.1
    )
    deg.i$layer <- sublayer.i
    deg.i$trans <- paste(sample.list[sample.i], sample.list[sample.i + 1], sep = ".")
    deg.i <- deg.i[deg.i$p_val_adj < 0.05,]
    deg.list.i <- c(deg.list.i, list(deg.i))
  }
  deg.df.i <- do.call("rbind", deg.list.i)
  deg.list <- c(deg.list, list(deg.df.i))
}
sublayer.stage.deg <- do.call("rbind", deg.list)

save(
  sublayer.stage.deg,
  file = "outputs/DEG/sublayer_stage_deg.RData"
)

# Sublayer-specific DEG ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
retina.st.domain$sublayer <- radial.anno[Cells(retina.st.domain), "sublayer"]

Idents(retina.st.domain) <- factor(retina.st.domain$sublayer, levels = radial.list)
sublayer.deg <- FindAllMarkers(
  retina.st.domain, test.use = "wilcox", only.pos = TRUE,
  logfc.threshold = 0.25, min.pct = 0.1
)
sublayer.deg <- sublayer.deg[sublayer.deg$p_val_adj < 0.05,]

save(
  sublayer.deg,
  file = "outputs/DEG/sublayer_deg.RData"
)

# Sublayer-specific DEG GO ####
deg.list <- unique(sublayer.deg$gene)
group <- data.frame(
  gene = sublayer.deg$gene,
  group = sublayer.deg$cluster
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
  file = "outputs/DEG/sublayer_deg_go.RData"
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
domain.deg <- domain.deg[domain.deg$p_val_adj < 0.05,]

save(
  domain.deg,
  file = "outputs/DEG/domain_deg.RData"
)

# Domain-specific DEG GO ####
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0,]
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
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0,]
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
