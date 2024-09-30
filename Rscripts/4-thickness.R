# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/ImageProcessing/radial_anno.RData")

retina.st.filter$layer <- radial.anno[Cells(retina.st.filter), "layer"]

# Calculate radial thickness ####
cal.radial.thickness <- function(
  label.PE, label.B1, label.B2, target, diameter
) {
  ## Calculate distance
  dist.PE <- sqrt(
    ((label.PE$x - target$x)^2) + ((label.PE$y - target$y)^2)
  )
  min.PE <- which.min(dist.PE)
  dist.B1 <- sqrt(
    ((label.B1$x - target$x)^2) + ((label.B1$y - target$y)^2)
  )
  min.B1 <- which.min(dist.B1)
  dist.B2 <- sqrt(
    ((label.B2$x - target$x)^2) + ((label.B2$y - target$y)^2)
  )
  min.B2 <- which.min(dist.B2)

  if (target$layer == "NBL") {
    dl <- sqrt(
      ((label.PE$x[min.PE] - label.B1$x[min.B1])^2) +
        ((label.PE$y[min.PE] - label.B1$y[min.B1])^2)
    )
  } else {
    dl <- sqrt(
      ((label.B1$x[min.B1] - label.B2$x[min.B2])^2) +
        ((label.B1$y[min.B1] - label.B2$y[min.B2])^2)
    )
  }
  dn <- sqrt(
    ((label.PE$x[min.PE] - label.B2$x[min.B2])^2) +
      ((label.PE$y[min.PE] - label.B2$y[min.B2])^2)
  )

  dl.thickness <- dl / diameter * 55
  dn.thickness <- dn / diameter * 55
  list(dl.thickness, dn.thickness)
}

# Load label & calculation ####
thickness.anno <- retina.st.filter@meta.data[, c("sample.name", "layer")]
thickness.anno$thickness.layer <- -1.0
thickness.anno$thickness.retina <- -1.0

mesh.list <- list()

for (sample.i in sample.list) {
  ## Load label
  label.i <- jsonlite::read_json(
    paste0("outputs/ImageProcessing/Annotation/", sample.i, ".json")
  )$shapes

  label.i.list <- list()
  for (point.i in seq_along(label.i)) {
    label.i.list[[point.i]] <- c(
      label.i[[point.i]]$label,
      label.i[[point.i]]$points[[1]][[1]],
      label.i[[point.i]]$points[[1]][[2]]
    )
  }
  label.i.df <- as.data.frame(t(as.data.frame(label.i.list)))
  colnames(label.i.df) <- c("label", "x", "y")
  label.i.df$x <- as.numeric(label.i.df$x)
  label.i.df$y <- as.numeric(label.i.df$y)

  label.i.mesh <- create.all.mesh(label.i.df)
  mesh.list[[sample.i]] <- label.i.mesh

  ## Load target
  idx.i <- retina.st.filter$sample.name == sample.i
  slice.i <- retina.st.filter$slice.name[idx.i][1]
  coo.i <- (
    retina.st.filter
      @images[[paste0("slice_", slice.i)]]
      @coordinates[, 4:5]
  )
  coo.i <- coo.i[Cells(retina.st.filter)[idx.i],]
  colnames(coo.i) <- c("y", "x")
  coo.i$layer <- retina.st.filter@meta.data[rownames(coo.i), "layer"]
  diameter.i <- (
    retina.st.filter
      @images[[paste0("slice_", slice.i)]]
      @scale.factors$spot
  )

  for (point.i in seq_len(nrow(coo.i))) {
    res.i <- cal.radial.thickness(
      label.i.mesh$PE, label.i.mesh$B1, label.i.mesh$B2,
      coo.i[point.i,], diameter.i
    )
    thickness.anno[rownames(coo.i)[point.i], "thickness.layer"] <- res.i[[1]]
    thickness.anno[rownames(coo.i)[point.i], "thickness.retina"] <- res.i[[2]]
  }
}

# Save data ####
save(
  thickness.anno,
  file = "outputs/ImageProcessing/thickness_anno.RData"
)

save(
  mesh.list,
  file = "outputs/ImageProcessing/mesh_for_layer.RData"
)