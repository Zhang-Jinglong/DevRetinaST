# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggforce)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate.RData")
load("outputs/Disease/kegg_go_gene.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

# Spatial plot function ####
spatial.gene.set.plot <- function(gene.set) {
  plot.data <- data.frame(
    value = colMeans(retina.st.domain@assays$Spatial@data[gene.set, ]),
    x = (
      retina.st.domain$coord.x +
        coord.add.1[retina.st.domain$sample.name, "x"]
    ),
    y = (
      retina.st.domain$coord.y +
        coord.add.1[retina.st.domain$sample.name, "y"]
    )
  )
  plot.data$value <- scale(plot.data$value)
  plot.data$value[plot.data$value > 3] <- 3
  plot.data$value[plot.data$value < -3] <- -3
  
  pic <- ggplot()
  for (sample.i in sample.list) {
    pic <- create.empty.bk(
      pic, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
    )
  }
  
  pic <- pic +
    geom_point(
      data = plot.data, size = 0.4,
      mapping = aes(x = x, y = y, color = value)
    ) +
    scale_color_gradient2(low = "#5D94A4", mid = "#FFFFFF", high = "#DA3B46") +
    theme_bw() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    coord_fixed()
  pic
}

# FAM spatial expression (Fig. 10a) ####
spatial.gene.set.plot(key.gene.set$FAM)

ggsave(
  "outputs/Visualization/fig-S10a.pdf",
  width = 10, height = 4
)

# FSVB spatial expression (Fig. 10b) ####
spatial.gene.set.plot(key.gene.set$FSVB)

ggsave(
  "outputs/Visualization/fig-S10b.pdf",
  width = 10, height = 4
)

# FSVC spatial expression (Fig. 10c) ####
spatial.gene.set.plot(key.gene.set$FSVC)

ggsave(
  "outputs/Visualization/fig-S10c.pdf",
  width = 10, height = 4
)

# VDM spatial expression (Fig. 10d) ####
spatial.gene.set.plot(key.gene.set$VDM)

ggsave(
  "outputs/Visualization/fig-S10d.pdf",
  width = 10, height = 4
)

# VEM spatial expression (Fig. 10e) ####
spatial.gene.set.plot(key.gene.set$VEM)

ggsave(
  "outputs/Visualization/fig-S10e.pdf",
  width = 10, height = 4
)

# VKM spatial expression (Fig. 10f) ####
spatial.gene.set.plot(key.gene.set$VKM)

ggsave(
  "outputs/Visualization/fig-S10f.pdf",
  width = 10, height = 4
)
