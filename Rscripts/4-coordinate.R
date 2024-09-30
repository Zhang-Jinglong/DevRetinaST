# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")

# Create new coordinate ####
spatial.coord.new <- data.frame(
  layer = radial.anno[Cells(retina.st.filter), "layer"],
  quantile = radial.anno[Cells(retina.st.filter), "quantile"],
  degree = degree.df[Cells(retina.st.filter), "degree"],
  row.names = Cells(retina.st.filter)
)
nbl.idx <- spatial.coord.new$layer == "NBL"
spatial.coord.new$quantile[nbl.idx] <- 7 - spatial.coord.new$quantile[nbl.idx]
gcl.idx <- spatial.coord.new$layer == "GCL"
spatial.coord.new$quantile[gcl.idx] <- 6 - spatial.coord.new$quantile[gcl.idx]

spatial.coord.new$coord.x <- (
  -spatial.coord.new$quantile * cos(spatial.coord.new$degree / 180 * pi)
)
spatial.coord.new$coord.y <- (
  spatial.coord.new$quantile * sin(spatial.coord.new$degree / 180 * pi)
)
spatial.coord.new <- spatial.coord.new[, 4:5]

# Save data ####
save(
  spatial.coord.new,
  file = "outputs/ImageProcessing/coordinate.RData"
)

# Extract H&E image & coordinate ####
library(ggplot2)
library(grid)
library(imager)
library(png)
library(jsonlite)

adjust.coo.img <- function(coo, img, scale) {
  img.new <- as.array(imresize(as.cimg(img), scale = scale))[, , 1,]

  scale.x <- dim(img.new)[2] / dim(img)[2]
  scale.y <- dim(img.new)[1] / dim(img)[1]
  coo.new <- coo[, c("coord.x", "coord.y")]
  coo.new$coord.x <- coo.new$coord.x * scale.x
  coo.new$coord.y <- dim(img.new)[1] - ((dim(img)[1] - coo.new$coord.y) * scale.y)

  coo.min <- as.integer(apply(coo.new, 2, min) - 10)
  coo.max <- as.integer(apply(coo.new, 2, max) + 10)
  coo.new <- t(t(coo.new) - coo.min)
  img.new <- img.new[(dim(img.new)[1] - coo.max[2]):(dim(img.new)[1] - coo.min[2]), coo.min[1]:coo.max[1],]

  return(list(coo.new, img.new))
}

## Own data ####
load("outputs/QualityControl/retina_st_filter.RData")

coo.list <- list()
img.list <- list()

img.obj.i <- retina.st.filter@images$slice_PCW9_12
img.i <- img.obj.i@image
img.i <- img.i[dim(img.i)[1]:1, ,]
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW9"],]
coo.i$coord.x <- coo.i$imagecol
coo.i$coord.y <- coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW9"]] <- as.data.frame(res.i[[1]])
img.list[["PCW9"]] <- res.i[[2]]

img.obj.i <- retina.st.filter@images$slice_PCW9_12
img.i <- img.obj.i@image
img.i <- img.i[, dim(img.i)[2]:1,]
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW12_early"],]
coo.i$coord.x <- dim(img.i)[2] - coo.i$imagecol
coo.i$coord.y <- dim(img.i)[1] - coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW12_early"]] <- as.data.frame(res.i[[1]])
img.list[["PCW12_early"]] <- res.i[[2]]

img.obj.i <- retina.st.filter@images$slice_PCW12
img.i <- img.obj.i@image
img.i <- img.i[, dim(img.i)[2]:1,]
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW12_late"],]
coo.i$coord.x <- dim(img.i)[2] - coo.i$imagecol
coo.i$coord.y <- dim(img.i)[1] - coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW12_late"]] <- as.data.frame(res.i[[1]])
img.list[["PCW12_late"]] <- res.i[[2]]

img.obj.i <- retina.st.filter@images$slice_PCW14
img.i <- img.obj.i@image
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW14"],]
coo.i$coord.x <- coo.i$imagecol
coo.i$coord.y <- dim(img.i)[1] - coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW14"]] <- as.data.frame(res.i[[1]])
img.list[["PCW14"]] <- res.i[[2]]

img.obj.i <- retina.st.filter@images$slice_PCW16
img.i <- img.obj.i@image
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW16"],]
coo.i$coord.x <- coo.i$imagecol
coo.i$coord.y <- dim(img.i)[1] - coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW16"]] <- as.data.frame(res.i[[1]])
img.list[["PCW16"]] <- res.i[[2]]

img.obj.i <- retina.st.filter@images$slice_PCW17
img.i <- img.obj.i@image
coo.i <- img.obj.i@coordinates[, c("imagerow", "imagecol")] * img.obj.i@scale.factors$lowres
coo.i <- coo.i[Cells(retina.st.filter)[retina.st.filter$sample.name == "PCW17"],]
coo.i$coord.x <- coo.i$imagecol
coo.i$coord.y <- dim(img.i)[1] - coo.i$imagerow
res.i <- adjust.coo.img(coo.i, img.i, 4 / (img.obj.i@scale.factors$spot * img.obj.i@scale.factors$lowres))
coo.list[["PCW17"]] <- as.data.frame(res.i[[1]])
img.list[["PCW17"]] <- res.i[[2]]

## GSE234035 ####
load("outputs/QualityControl/retina_st_pub.RData")

pub.coo.list <- list()
pub.img.list <- list()

for (prefix.i in unique(retina.st.pub$image.prefix)) {
  coo.i <- retina.st.pub@meta.data[retina.st.pub$image.prefix == prefix.i, c("x", "y")]
  img.i <- readPNG(paste0(
    "data/ST/GSE234035/GSE234035_RAW/", prefix.i, "_tissue_lowres_image.png"
  ))
  coo.i$coord.x <- coo.i$x
  coo.i$coord.y <- dim(img.i)[1] - coo.i$y

  if (prefix.i == "GSM7441454_14749_13PCW_A") {
    img.i <- img.i[dim(img.i)[1]:1, dim(img.i)[2]:1,]
    coo.i$coord.x <- dim(img.i)[2] - coo.i$coord.x
    coo.i$coord.y <- dim(img.i)[1] - coo.i$coord.y
  }

  if (strsplit(prefix.i, "_")[[1]][3] != "8PCW") {
    img.i <- aperm(img.i, perm = c(2, 1, 3))
    coo.i.copy <- coo.i
    coo.i$coord.x <- dim(img.i)[1] - coo.i.copy$coord.y
    coo.i$coord.y <- dim(img.i)[2] - coo.i.copy$coord.x
  } else {
    img.i <- img.i[dim(img.i)[1]:1, dim(img.i)[2]:1,]
    coo.i$coord.x <- dim(img.i)[2] - coo.i$coord.x
    coo.i$coord.y <- dim(img.i)[1] - coo.i$coord.y
  }

  res.i <- adjust.coo.img(coo.i, img.i, 11 / 13)
  pub.coo.list[[prefix.i]] <- as.data.frame(res.i[[1]])
  pub.img.list[[prefix.i]] <- res.i[[2]]
}

## Save data ####
save(
  coo.list, img.list, pub.coo.list, pub.img.list,
  file = "outputs/ImageProcessing/coordinate_HE.RData"
)
