# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# CellPhoneDB ####
## Convert input data ####
file.symlink(
  from = paste0(
    getwd(), "/outputs/Deconvolution/Cell2location/retina_sc.csv"
  ),
  to = "outputs/CCC/CellPhoneDB/retina_sc.csv"
)
file.symlink(
  from = paste0(
    getwd(), "/outputs/Deconvolution/Cell2location/retina_sc_meta.csv"
  ),
  to = "outputs/CCC/CellPhoneDB/retina_sc_meta.csv"
)

## python 9-cellphonedb.ipynb ####
## outputs: "outputs/CCC/CellPhoneDB/"

## Matching with COMMOT database ####
### Build database mapping ####
cpdb.db.filter <- read.csv(
  "outputs/CCC/COMMOT/retina_db.csv", row.names = 1
)
rownames(cpdb.db.filter) <- paste0("commot_", 1:nrow(cpdb.db.filter))
cpdb.db.filter.list <- list()
for (db.i in rownames(cpdb.db.filter)) {
  db.gene.i <- c(
    unlist(strsplit(cpdb.db.filter[db.i, "X0"], split = "_")),
    unlist(strsplit(cpdb.db.filter[db.i, "X1"], split = "_"))
  )
  cpdb.db.filter.list[[db.i]] <- db.gene.i
}

cpdb.mapping <- read.csv(
  "outputs/CCC/CellPhoneDB/statistical_analysis_deconvoluted.txt",
  sep = "\t"
)[, 1:6]
cpdb.mapping.commot <- data.frame(
  cpi.id = unique(cpdb.mapping$id_cp_interaction),
  commot = NA
)
rownames(cpdb.mapping.commot) <- cpdb.mapping.commot$cpi.id
for (cpi.i in rownames(cpdb.mapping.commot)) {
  cpi.idx.i <- cpdb.mapping$id_cp_interaction == cpi.i
  cpi.gene.i <- unique(cpdb.mapping$gene_name[cpi.idx.i])
  
  for (db.i in names(cpdb.db.filter.list)) {
    db.gene.i <- cpdb.db.filter.list[[db.i]]
    
    if (
      (length(cpi.gene.i) == length(db.gene.i)) &
      (length(intersect(cpi.gene.i, db.gene.i)) == length(db.gene.i))
    ) {
      cpdb.mapping.commot[cpi.i, "commot"] <- db.i
      cpdb.db.filter.list[[db.i]] <- NULL
      break
    }
  }
}
cpdb.mapping.commot.idx <- !is.na(cpdb.mapping.commot$commot)
cpdb.mapping.commot <- cpdb.mapping.commot[cpdb.mapping.commot.idx, ]
rownames(cpdb.mapping.commot) <- cpdb.mapping.commot$commot

cpdb.db.filter <- cpdb.db.filter[cpdb.mapping.commot$commot, ]
rm(cpdb.db.filter.list, cpdb.mapping.commot.idx)

### Matching result ####
#### Create gene name reference ####
cpdb.mapping.list <- c()
single.list <- unique(
  cpdb.mapping[cpdb.mapping$is_complex == "False", 1:2]
)
for (single.i in 1:nrow(single.list)) {
  single.gene.i <- single.list$gene_name[single.i]
  cpdb.mapping.list[single.list$uniprot[single.i]] <- single.gene.i
}
complex.list <- unique(
  cpdb.mapping$complex_name[cpdb.mapping$is_complex == "True"]
)
for (complex.i in complex.list) {
  complex.gene.i <- unique(
    cpdb.mapping$gene_name[cpdb.mapping$complex_name == complex.i]
  )
  if (length(complex.gene.i) > 1) {
    complex.gene.i <- paste0(complex.gene.i, collapse = "_")
  }
  
  cpdb.mapping.list[complex.i] <- complex.gene.i
}
rm(single.list, complex.list, cpdb.mapping)

#### Load results ####
cpdb.means <- read.csv(
  "outputs/CCC/CellPhoneDB/statistical_analysis_means.txt",
  sep = "\t", row.names = 1
)[cpdb.mapping.commot$cpi.id, c(-1, -6, -10)]
cpdb.p.val <- read.csv(
  "outputs/CCC/CellPhoneDB/statistical_analysis_pvalues.txt",
  sep = "\t", row.names = 1
)[cpdb.mapping.commot$cpi.id, c(-1, -6, -10)]

#### Transfer gene name ####
cpdb.means.split <- sapply(strsplit(cpdb.means$partner_a, split = ":"), "[", -1)
cpdb.means.a <- c()
for (cpdb.means.i in cpdb.means.split) {
  cpdb.means.a <- c(cpdb.means.a, paste(cpdb.means.i, collapse = ":"))
}
cpdb.means.a <- cpdb.mapping.list[cpdb.means.a]
cpdb.means.split <- sapply(strsplit(cpdb.means$partner_b, split = ":"), "[", -1)
cpdb.means.b <- c()
for (cpdb.means.i in cpdb.means.split) {
  cpdb.means.b <- c(cpdb.means.b, paste(cpdb.means.i, collapse = ":"))
}
cpdb.means.b <- cpdb.mapping.list[cpdb.means.b]
rm(cpdb.means.split)

cpdb.means$gene_a <- cpdb.means.a
cpdb.means$gene_b <- cpdb.means.b
cpdb.means$annotation_strategy <- cpdb.mapping.commot$anno
cpdb.p.val$gene_a <- cpdb.means.a
cpdb.p.val$gene_b <- cpdb.means.b
cpdb.p.val$annotation_strategy <- cpdb.mapping.commot$anno

#### Filtering L-R pairs ####
cpdb.ava.idx <- (
  ((cpdb.means$receptor_a == "True") & (cpdb.means$receptor_b == "True")) |
    ((cpdb.means$receptor_a != "True") & (cpdb.means$receptor_b != "True"))
)
cpdb.means <- cpdb.means[!cpdb.ava.idx, ]
cpdb.p.val <- cpdb.p.val[!cpdb.ava.idx, ]
cpdb.gene.a <- cpdb.means$gene_a
cpdb.gene.b <- cpdb.means$gene_b
swap.idx <- cpdb.means$receptor_a == "True"
tmp <- cpdb.gene.a[swap.idx]
cpdb.gene.a[swap.idx] <- cpdb.gene.b[swap.idx]
cpdb.gene.b[swap.idx] <- tmp
cpdb.means$gene_a <- cpdb.gene.a
cpdb.means$gene_b <- cpdb.gene.b
cpdb.p.val$gene_a <- cpdb.gene.a
cpdb.p.val$gene_b <- cpdb.gene.b

cpdb.dup.idx <- duplicated(cpdb.means[, c("gene_a", "gene_b")])
cpdb.means <- cpdb.means[!cpdb.dup.idx, ]
cpdb.p.val <- cpdb.p.val[!cpdb.dup.idx, ]

cpdb.ava.idx <- rowSums(cpdb.p.val[, 7:106] < 0.05) > 0
cpdb.means.filter <- cpdb.means[cpdb.ava.idx, ]
cpdb.p.val.filter <- cpdb.p.val[cpdb.ava.idx, ]

write.csv(
  cpdb.means.filter[, c(-1, -2, -5, -6)],
  file = "outputs/CCC/CellPhoneDB/cpdb_filter_means.csv"
)
write.csv(
  cpdb.p.val.filter[, c(-1, -2, -5, -6)],
  file = "outputs/CCC/CellPhoneDB/cpdb_filter_pvalues.csv"
)

### Save filtered database for COMMOT ####
cpi.idx <- cpdb.mapping.commot$cpi.id %in% rownames(cpdb.means)
cpdb.mapping.commot <- cpdb.mapping.commot[cpi.idx, ]
cpdb.mapping.commot$gene.a <- cpdb.means[cpdb.mapping.commot$cpi.id, "gene_a"]
cpdb.mapping.commot$gene.b <- cpdb.means[cpdb.mapping.commot$cpi.id, "gene_b"]

cpdb.db.filter <- cpdb.db.filter[cpdb.mapping.commot$commot, ]
cpdb.db.filter$X0 <- cpdb.mapping.commot$gene.a
cpdb.db.filter$X1 <- cpdb.mapping.commot$gene.b

cpdb.db.filter$X2 <- ""
rownames(cpdb.db.filter) <- seq(1:nrow(cpdb.db.filter)) - 1

write.csv(
  cpdb.db.filter, quote = FALSE,
  file = "outputs/CCC/COMMOT/retina_db_filter.csv"
)

# COMMOT ####
## Convert input data ####
st.count <- data.frame(
  retina.st.domain@assays$Spatial@counts,
  check.names = FALSE
)
write.csv(
  st.count,
  file = "outputs/CCC/COMMOT/retina_st.csv"
)

st.metadata <- data.frame(
  sample = retina.st.domain$sample.name,
  domain = retina.st.domain$domain,
  coord.x = retina.st.domain$coord.x,
  coord.y = retina.st.domain$coord.y,
  row.names = Cells(retina.st.domain)
)
write.csv(
  st.metadata,
  file = "outputs/CCC/COMMOT/retina_st_meta.csv"
)

## python 9-commot.py ####
## outputs: "outputs/CCC/COMMOT/"
