# Packages & settings ####
library(Seurat)
library(GeneOverlap)
library(VennDiagram)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate.RData")
load("outputs/DEG/domain_deg.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Load TF/RBP list ####
act.gene <- rownames(retina.st.domain)
act.gene <- act.gene[
  rowSums(retina.st.domain@assays$Spatial@counts > 0) >
    (ncol(retina.st.domain) * 0.1)
]

tf.info <- readxl::read_xlsx(
  "outputs/Reference/Cell_2018_TF_List/mmc2.xlsx",
  sheet = 2, range = "B2:D2767"
)
tf.list <- tf.info$Name[tf.info$...3 == "Yes"]
tf.list <- intersect(tf.list, act.gene)

rbp.info <- readxl::read_xlsx(
  paste0(
    "outputs/Reference/NatRevMolCellBiol_2018_RBP_List/",
    "41580_2018_BFnrm2017130_MOESM121_ESM.xlsx"
  ), sheet = "Hs", range = "C2:C1395"
)
rbp.list <- rbp.info$...1
rbp.list <- intersect(rbp.list, act.gene)

save(
  tf.list, rbp.list,
  file = "outputs/Development/reg_gene_ava.RData"
)

# Key TF/RBP ####
domain.sub.list <- domain.list[1:8]
trans.list <- data.frame(
  source = sample.list[1:5],
  target = sample.list[2:6]
)

deg.list <- list()
for (domain.i in domain.sub.list) {
  retina.st.sub <- subset(
    retina.st.domain,
    subset = domain == domain.i
  )
  Idents(retina.st.sub) <- retina.st.sub$sample.name
  
  for (trans.i in 1:nrow(trans.list)) {
    source.i <- trans.list[trans.i, "source"]
    target.i <- trans.list[trans.i, "target"]
    
    deg.i <- FindMarkers(
      retina.st.sub, test.use = "wilcox",
      ident.1 = target.i, ident.2 = source.i,
      features = unique(c(tf.list, rbp.list)),
      min.pct = 0.1, logfc.threshold = 0.25,
      only.pos = TRUE, verbose = FALSE
    )
    deg.i$p_val_adj <- p.adjust(deg.i$p_val, method = "bonferroni")
    
    if (nrow(deg.i) > 0) {
      deg.i$gene <- rownames(deg.i)
      deg.i$domain <- domain.i
      deg.i$source <- source.i
      deg.i$target <- target.i
      deg.list <- c(deg.list, list(deg.i))
    }
  }
}
reg.gene.df <- do.call("rbind", deg.list)
reg.gene.df <- reg.gene.df[reg.gene.df$p_val_adj < 0.05, ]

both.list <- intersect(tf.list, rbp.list)
tf.uni.list <- setdiff(tf.list, both.list)
rbp.uni.list <- setdiff(rbp.list, both.list)
reg.gene.df$type <- "Both"
reg.gene.df$type[reg.gene.df$gene %in% tf.uni.list] <- "TF only"
reg.gene.df$type[reg.gene.df$gene %in% rbp.uni.list] <- "RBP only"

save(
  reg.gene.df,
  file = "outputs/Development/reg_gene_df.RData"
)
