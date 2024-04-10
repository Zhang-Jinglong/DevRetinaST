# Packages & settings ####
library(Seurat)
library(SPARK)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")

# SPARK-X ####
svg.list <- list()
for (sample.i in sample.list) {
  obj.i <- subset(retina.st.filter, subset = sample.name == sample.i)
  img.i <- obj.i@images[[paste0("slice_", obj.i$slice.name[1])]]
  
  sp.count.i <- obj.i@assays$Spatial@counts
  gene.ava <- rowSums(sp.count.i > 0) >= 10
  sp.count.i <- sp.count.i[gene.ava, ]
  
  info.i <- cbind.data.frame(
    x = img.i@coordinates$imagerow,
    y = img.i@coordinates$imagecol
  )
  rownames(info.i) <- colnames(sp.count.i)
  location.i <- as.matrix(info.i)
  
  spark.x.i <- sparkx(
    sp.count.i, location.i, numCores = 4, option = "mixture"
  )
  spark.x.i$res_mtest$adjustedPval <- p.adjust(
    spark.x.i$res_mtest$combinedPval, method = "bonferroni"
  )
  
  idx.i <- spark.x.i$res_mtest$adjustedPval < 0.05
  svg.list[[sample.i]] <- data.frame(
    gene.symbol = rownames(spark.x.i$res_mtest)[idx.i],
    p.adj = spark.x.i$res_mtest$adjustedPval[idx.i]
  )
}

# Save data ####
save(svg.list, file = "outputs/SVG/svg_list.RData")
