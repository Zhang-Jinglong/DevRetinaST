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

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

# VEGFA spatial expression (Fig. S9a) ####
plot.data <- data.frame(
  value = retina.st.domain@assays$Spatial@data["VEGFA", ],
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

pic.a <- ggplot()
for (sample.i in sample.list) {
  pic.a <- create.empty.bk(
    pic.a, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.a +
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

ggsave(
  "outputs/Visualization/fig-S9a.pdf",
  width = 10, height = 4
)

# PGF spatial expression (Fig. S9b) ####
plot.data <- data.frame(
  value = retina.st.domain@assays$Spatial@data["PGF", ],
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

pic.b <- ggplot()
for (sample.i in sample.list) {
  pic.b <- create.empty.bk(
    pic.b, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.b +
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

ggsave(
  "outputs/Visualization/fig-S9b.pdf",
  width = 10, height = 4
)
