# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
radial.anno <- data.frame(
  layer = rep("", ncol(retina.st.filter)),
  quantile = rep(Inf, ncol(retina.st.filter)),
  row.names = Cells(retina.st.filter)
)

# Calculate radial annotation ####
cal.radial.anno <- function(
    label.PE, label.B1, label.B2, target
) {
  ## Calculate distance
  dist.PE <- sqrt(
    ((label.PE$x - target$x) ^ 2) + ((label.PE$y - target$y) ^ 2)
  )
  min.PE <- which.min(dist.PE)
  dp <- dist.PE[min.PE]
  
  dist.B1 <- sqrt(
    ((label.B1$x - target$x) ^ 2) + ((label.B1$y - target$y) ^ 2)
  )
  min.B1 <- which.min(dist.B1)
  d1 <- dist.B1[min.B1]
  
  dist.B2 <- sqrt(
    ((label.B2$x - target$x) ^ 2) + ((label.B2$y - target$y) ^ 2)
  )
  min.B2 <- which.min(dist.B2)
  d2 <- dist.B2[min.B2]
  
  da.p2 <- sqrt(
    ((label.PE$x[min.PE] - label.B2$x[min.B2]) ^ 2) +
      ((label.PE$y[min.PE] - label.B2$y[min.B2]) ^ 2)
  )
  da.12 <- sqrt(
    ((label.B1$x[min.B1] - label.B2$x[min.B2]) ^ 2) +
      ((label.B1$y[min.B1] - label.B2$y[min.B2]) ^ 2)
  )
  
  ## Annotation
  if ((d2 > da.p2) & (d2 > dp)) {
    target.layer <- "NBL"
    target.qt <- - dp / da.p2
  } else if ((dp > da.p2) & (dp > d2 )) {
    target.layer <- "GCL"
    target.qt <- 1 + (d2 / da.p2)
  } else if (d2 > da.12) {
    target.layer <- "NBL"
    target.qt <- dp / (dp + d1)
  } else {
    target.layer <- "GCL"
    target.qt <- d1 / (d1 + d2)
  }
  
  data.frame(layer = target.layer, quantile = target.qt)
}

# Load label & calculation ####
for (sample.i in sample.list) {
  ## Load label
  label.i <- jsonlite::read_json(
    paste0("outputs/ImageProcessing/Annotation/", sample.i, ".json")
  )$shapes
  
  label.i.list <- list()
  for (point.i in 1:length(label.i)) {
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
  
  ## Load target
  idx.i <- retina.st.filter$sample.name == sample.i
  slice.i <- retina.st.filter$slice.name[idx.i][1]
  coo.i <- (
    retina.st.filter
    @images[[paste0("slice_", slice.i)]]
    @coordinates[, 4:5]
  )
  coo.i <- coo.i[Cells(retina.st.filter)[idx.i], ]
  colnames(coo.i) <- c("y", "x")
  
  for (point.i in 1:nrow(coo.i)) {
    res.i <- cal.radial.anno(
      label.i.mesh$PE, label.i.mesh$B1, label.i.mesh$B2,
      coo.i[point.i, ]
    )
    radial.anno[rownames(coo.i)[point.i], "layer"] <- res.i$layer
    radial.anno[rownames(coo.i)[point.i], "quantile"] <- res.i$quantile
  }
}

# Quantile normalization ####
radial.nbl <- radial.anno[radial.anno$layer == "NBL", ]
q99.nbl <- quantile(radial.nbl$quantile, 1 - 0.99)
q99.nbl.idx <- radial.nbl$quantile < q99.nbl
radial.nbl$quantile[q99.nbl.idx] <- q99.nbl
radial.nbl$quantile <- rank(radial.nbl$quantile) / nrow(radial.nbl)

radial.gcl <- radial.anno[radial.anno$layer == "GCL", ]
q99.gcl <- quantile(radial.gcl$quantile, 0.99)
q99.gcl.idx <- radial.gcl$quantile > q99.gcl
radial.gcl$quantile[q99.gcl.idx] <- q99.gcl
radial.gcl$quantile <- rank(radial.gcl$quantile) / nrow(radial.gcl)

radial.anno <- rbind(radial.nbl, radial.gcl)

# Sub-layer annotation ####
sub.idx.1 <- (radial.anno$layer == "GCL") &
  (radial.anno$quantile < 0.5)
radial.anno$sublayer[sub.idx.1] <- "GCL 1/2"

sub.idx.2 <- (radial.anno$layer == "GCL") &
  (radial.anno$quantile >= 0.5)
radial.anno$sublayer[sub.idx.2] <- "GCL 2/2"

sub.idx.3 <- (radial.anno$layer == "NBL") &
  (radial.anno$quantile < 0.25)
radial.anno$sublayer[sub.idx.3] <- "NBL 1/4"

sub.idx.4 <- (radial.anno$layer == "NBL") &
  (radial.anno$quantile < 0.5) &
  (radial.anno$quantile >= 0.25)
radial.anno$sublayer[sub.idx.4] <- "NBL 2/4"

sub.idx.5 <- (radial.anno$layer == "NBL") &
  (radial.anno$quantile < 0.75) &
  (radial.anno$quantile >= 0.5)
radial.anno$sublayer[sub.idx.5] <- "NBL 3/4"

sub.idx.6 <- (radial.anno$layer == "NBL") &
  (radial.anno$quantile >= 0.75)
radial.anno$sublayer[sub.idx.6] <- "NBL 4/4"

# Save data ####
save(
  radial.anno,
  file = "outputs/ImageProcessing/radial_anno.RData"
)
