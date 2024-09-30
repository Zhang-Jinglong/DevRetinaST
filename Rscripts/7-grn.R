# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/SVG/svg_list.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)
gene.list <- unique(do.call("rbind", svg.list)$gene.symbol)

mat <- as.data.frame(retina.st.domain@assays$Spatial@counts[gene.list,])
write.csv(mat, file = "outputs/Development/SCRIPro/data.csv")
meta <- data.frame(
  spot = colnames(retina.st.domain),
  sample = retina.st.domain$sample.name,
  domain = retina.st.domain$domain
)
write.csv(meta, file = "outputs/Development/SCRIPro/meta.csv")

# Load TF ####
load("outputs/Development/reg_gene_df.RData")

tf.list <- unique(reg.gene.df$gene[reg.gene.df$type != "RBP only"])

write.csv(tf.list, file = "outputs/Development/SCRIPro/reg_tf_list.csv")
