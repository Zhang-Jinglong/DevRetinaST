# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")

# Load disease gene set ####
retnet.table <- readxl::read_excel(
  "outputs/Reference/RetNet_20221007.xlsx",
  sheet = "Merge"
)
retnet <- list()
for (dis.i in 1:nrow(retnet.table)) {
  gene.i <- unique(
    strsplit(retnet.table$Genes[dis.i], split = ",")[[1]]
  )
  retnet[[retnet.table$Disease[dis.i]]] <- intersect(
    gene.i, rownames(retina.st.domain)
  )
}

save(retnet, file = "outputs/Disease/retnet.RData")

# KEGG & GO gene sets ####
kegg.gene.set <- read.gmt(
  "outputs/Reference/c2.cp.kegg.v2023.1.Hs.symbols.gmt"
)
kegg.idx <- kegg.gene.set$gene %in% rownames(retina.st.domain)
kegg.gene.set <- kegg.gene.set[kegg.idx, ]
go.gene.set <- read.gmt(
  "outputs/Reference/c5.go.bp.v2023.1.Hs.symbols.gmt"
)
go.idx <- go.gene.set$gene %in% rownames(retina.st.domain)
go.gene.set <- go.gene.set[go.idx, ]

key.gene.set <- list(
  FAM = kegg.gene.set[
    kegg.gene.set$term == "KEGG_FATTY_ACID_METABOLISM", "gene"
  ],
  FSVB = go.gene.set[
    go.gene.set$term == "GOBP_FAT_SOLUBLE_VITAMIN_BIOSYNTHETIC_PROCESS", "gene"
  ],
  FSVC = go.gene.set[
    go.gene.set$term == "GOBP_FAT_SOLUBLE_VITAMIN_CATABOLIC_PROCESS", "gene"
  ],
  VAM = go.gene.set[
    go.gene.set$term == "GOBP_VITAMIN_A_METABOLIC_PROCESS", "gene"
  ],
  VDM = go.gene.set[
    go.gene.set$term == "GOBP_VITAMIN_D_METABOLIC_PROCESS", "gene"
  ],
  VEM = go.gene.set[
    go.gene.set$term == "GOBP_VITAMIN_E_METABOLIC_PROCESS", "gene"
  ],
  VKM = go.gene.set[
    go.gene.set$term == "GOBP_VITAMIN_K_METABOLIC_PROCESS", "gene"
  ]
)

save(key.gene.set, file = "outputs/Disease/kegg_go_gene.RData")
