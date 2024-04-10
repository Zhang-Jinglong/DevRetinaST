# Packages & settings ####
library(Seurat)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/MergeData/retina_st.RData")

for (slice.i in slice.list) {
  retina.st[[slice.i]]$tangential.degree <- NA
  for (sample.i in unique(retina.st[[slice.i]]$sample.name)) {
    obj <- subset(
      retina.st[[slice.i]],
      subset = sample.name == sample.i
    )
    
    tangential <- read.csv(
      paste0("data/ST/", slice.i, "/outs/annotation/TangentialAnno.csv"),
      row.names = 1
    )
    obj$tangential.anno <- NA
    obj.names <- sapply(strsplit(Cells(obj), "[.]"), "[", 2)
    obj$tangential.anno <- tangential[obj.names, "TangentialAnno"]
    
    tmp.flag <- 1
    if (slice.i == "PCW9_12" & sample.i == "PCW9") {
      tmp.flag <- -1
    }
    
    # For superior area ####
    spot.superior.list <- Cells(obj)[obj$tangential.anno == "Superior"]
    spot.superior.center <- as.numeric(obj@images[[1]]@coordinates[
      obj$tangential.anno == "Superior-center", c("imagerow", "imagecol")
    ])
    spot.superior.list <- data.frame(
      x = spot.superior.center[1] - obj@images[[1]]@coordinates[
        spot.superior.list, "imagerow"
      ],
      y = spot.superior.center[2] - obj@images[[1]]@coordinates[
        spot.superior.list, "imagecol"
      ],
      row.names = spot.superior.list
    )
    spot.superior.list$degree <- atan(
      spot.superior.list$x / spot.superior.list$y
    ) * 180 / pi
    spot.superior.list$degree <- spot.superior.list$degree + (
      spot.superior.list$y < 0
    ) * 180 * tmp.flag
    
    spot.superior.list$degree <- spot.superior.list$degree * tmp.flag
    
    degree.min <- min(spot.superior.list[
      rownames(spot.superior.list), "degree"
    ])
    degree.max <- max(spot.superior.list[
      rownames(spot.superior.list), "degree"
    ])
    degree.cut <- degree.min
    spot.superior.list$degree <- (
      spot.superior.list$degree - degree.cut
    ) / (degree.max - degree.min) * (180 - 45)
    if (
      (slice.i == "PCW12" & sample.i == "PCW12_late") |
      (slice.i == "PCW9_12" & sample.i == "PCW12_early")
    ) {
      spot.superior.list$degree <- 180 - 45 * tmp.flag -
        spot.superior.list$degree
    }
    
    # For inferior area ####
    spot.inferior.list <- Cells(obj)[obj$tangential.anno == "Inferior"]
    spot.inferior.center <- as.numeric(obj@images[[1]]@coordinates[
      obj$tangential.anno == "Inferior-center", c("imagerow", "imagecol")
    ])
    spot.inferior.list <- data.frame(
      x = spot.inferior.center[1] - obj@images[[1]]@coordinates[
        spot.inferior.list, "imagerow"
      ],
      y = spot.inferior.center[2] - obj@images[[1]]@coordinates[
        spot.inferior.list, "imagecol"
      ],
      row.names = spot.inferior.list
    )
    spot.inferior.list$degree <- atan(
      spot.inferior.list$x / spot.inferior.list$y
    ) * 180 / pi
    spot.inferior.list$degree <- spot.inferior.list$degree - (
      spot.inferior.list$y < 0
    ) * 180 * tmp.flag
    
    spot.inferior.list$degree <- spot.inferior.list$degree * tmp.flag
    
    degree.min <- min(spot.inferior.list[
      rownames(spot.inferior.list), "degree"
    ])
    degree.max <- max(spot.inferior.list[
      rownames(spot.inferior.list), "degree"
    ])
    degree.cut <- degree.max
    spot.inferior.list$degree <- (
      spot.inferior.list$degree - degree.cut
    ) / (degree.max - degree.min) * (180 - 45)
    if (
      (slice.i == "PCW12" & sample.i == "PCW12_late") |
      (slice.i == "PCW9_12" & sample.i == "PCW12_early")
    ) {
      spot.inferior.list$degree <- - 180 + 45 * tmp.flag -
        spot.inferior.list$degree
    }
    
    # Degree annotation ####
    spot.tangential <- rbind(spot.superior.list, spot.inferior.list)
    retina.st[[slice.i]]$tangential.degree[
      rownames(spot.tangential)
    ] <- spot.tangential$degree
  }
}

tangential.degree.list <- list()
for (slice.i in slice.list) {
  degree.i <- retina.st[[slice.i]]$tangential.degree
  degree.i <- degree.i[!is.na(degree.i)]
  tangential.degree.list[[slice.i]] <- as.data.frame(degree.i)
}
degree.df <- do.call("rbind", tangential.degree.list)
new.names <- paste0(
  sapply(strsplit(rownames(degree.df), "[.]"), "[", 2), ".",
  sapply(strsplit(rownames(degree.df), "[.]"), "[", 3)
)
rownames(degree.df) <- new.names
colnames(degree.df) <- "degree"

# Save data ####
save(
  degree.df,
  file = "outputs/ImageProcessing/tangential_anno.RData"
)
