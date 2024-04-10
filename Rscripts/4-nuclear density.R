# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")

nuclear.density <- list()
for (slice.i in slice.list) {
  idx.i <- retina.st.filter$slice.name == slice.i
  nuclear.density.i <- data.frame(
    nuclear.density = rep(0.0, sum(idx.i)),
    row.names = Cells(retina.st.filter)[idx.i]
  )
  
  img.i <- jpeg::readJPEG(
    paste0("outputs/ImageProcessing/Binary/", slice.i, "_norm.jpg")
  )
  radius.i <- (
    retina.st.filter
    @images[[paste0("slice_", slice.i)]]
    @scale.factors$fiducial
  ) / 2
  coo.i <- (
    retina.st.filter
    @images[[paste0("slice_", slice.i)]]
    @coordinates[, 4:5]
  )
  coo.i <- coo.i[rownames(nuclear.density.i), ]
  
  for (spot.i in 1:nrow(coo.i)) {
    row.min <- as.integer(coo.i[spot.i, 1] - radius.i) - 1
    row.max <- as.integer(coo.i[spot.i, 1] + radius.i) + 1
    col.min <- as.integer(coo.i[spot.i, 2] - radius.i) - 1
    col.max <- as.integer(coo.i[spot.i, 2] + radius.i) + 1
    
    count.i <- 0
    for (row.i in row.min:row.max) {
      for (col.i in col.min:col.max) {
        if (img.i[row.i, col.i] < 0.5) {
          dist.i <- sqrt(
            ((row.i - coo.i[spot.i, 1]) ^ 2) +
              ((col.i - coo.i[spot.i, 2]) ^ 2)
          )
          if (dist.i < radius.i) { count.i = count.i + 1}
        }
      }
    }
    density.i <- count.i / pi / radius.i / radius.i
    
    name.i <- rownames(coo.i)[spot.i]
    nuclear.density.i[name.i, "nuclear.density"] <- density.i
  }
  
  nuclear.density[[slice.i]] <- nuclear.density.i
}

nuclear.density.df <- do.call("rbind", nuclear.density)
new.names <- paste0(
  sapply(strsplit(rownames(nuclear.density.df), "[.]"), "[", 2), ".",
  sapply(strsplit(rownames(nuclear.density.df), "[.]"), "[", 3)
)
rownames(nuclear.density.df) <- new.names

# Save data ####
save(
  nuclear.density.df,
  file = "outputs/ImageProcessing/nuclear_density.RData"
)
