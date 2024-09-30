# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggforce)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/Deconvolution/decon_results.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.filter)
retina.st.filter$layer <- radial.anno[spot.list, "layer"]
retina.st.filter$quantile <- radial.anno[spot.list, "quantile"]
retina.st.filter$tangential <- degree.df[spot.list, "degree"]
retina.st.filter$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.filter$coord.y <- spatial.coord.new[spot.list, "coord.y"]

# Spatial-temporal cell-type distribution (Fig. S5a) ####
coord.add.t <- data.frame(
  x = rep(0, 4),
  y = seq(0, 3, 1) * (-15),
  row.names = sort(cell.type.list[7:10])
)

pic.a <- ggplot()
for (type.i in sort(cell.type.list[7:10])) {
  plot.data.i <- data.frame(
    value = pred.cell2loction[type.i, Cells(retina.st.filter)],
    x = (
      retina.st.filter$coord.x +
        coord.add.t[type.i, "x"] +
        coord.add.1[retina.st.filter$sample.name, "x"]
    ),
    y = (
      retina.st.filter$coord.y +
        coord.add.t[type.i, "y"] +
        coord.add.1[retina.st.filter$sample.name, "y"]
    )
  )
  
  for (sample.i in sample.list) {
    pic.a <- create.empty.bk(
      pic.a,
      coord.add.t[type.i, "x"] + coord.add.1[sample.i, "x"],
      coord.add.t[type.i, "y"] + coord.add.1[sample.i, "y"]
    )
  }
  
  pic.a <- pic.a +
    geom_point(
      data = plot.data.i, size = 0.4,
      mapping = aes(x = x, y = y, color = value)
    ) +
    scale_color_gradient(low = "#FFFFFF", high = cell.type.colors[type.i]) +
    new_scale_color()
}
pic.a +
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
  "outputs/Visualization/fig-S5a.pdf",
  width = 10, height = 8
)

# Radial cell type distribution (Fig. S5b) ####
meta <- retina.st.filter@meta.data
meta$quantile[meta$layer == "GCL"] <- meta$quantile[meta$layer == "GCL"] + 1

meta.list <- list()
for (type.i in cell.type.list) {
  meta.list[[type.i]] <- data.frame(
    type = rep(type.i, nrow(meta)),
    num = correct.cell2loction[type.i, rownames(meta)],
    loc = meta$quantile,
    sample = meta$sample.name,
    row.names = rownames(meta)
  )
}
meta.df <- do.call("rbind", meta.list)

pic.list <- list()
for (sample.i in sample.list) {
  pic.list[[sample.i]] <- ggplot(meta.df[meta.df$sample == sample.i, ]) +
    geom_smooth(
      aes(x = loc, y = num, group = type, color = type),
      method = "glm", se = FALSE, formula = y ~ poly(x, 3),
      method.args = list(family = quasibinomial)
    ) + theme_classic() +
    scale_color_manual(values = cell.type.colors) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0),
      labels = c(radial.list, "Border")
    ) +
    scale_y_continuous(
      breaks = seq(0, 10) / 10,
      labels = seq(0, 10) * 10
    ) +
    coord_cartesian(ylim = c(0, 0.65)) +
    geom_vline(xintercept = c(0.25, 0.5, 0.75, 1.0, 1.5), linetype = 3) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 0.1),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

pic.list[["PCW9"]] | pic.list[["PCW12_early"]] | pic.list[["PCW12_late"]] |
  pic.list[["PCW14"]] | pic.list[["PCW16"]] | pic.list[["PCW17"]]

ggsave(
  filename = "outputs/Visualization/fig-S5b.pdf",
  width = 20, height = 4
)

# Tangential cell type distribution (Fig. S5c) ####
meta <- retina.st.filter@meta.data
meta$tangential <- abs(meta$tangential)

meta.list <- list()
for (type.i in cell.type.list) {
  meta.list[[type.i]] <- data.frame(
    type = rep(type.i, nrow(meta)),
    num = correct.cell2loction[type.i, rownames(meta)],
    loc = meta$tangential,
    sample = meta$sample.name,
    row.names = rownames(meta)
  )
}
meta.df <- do.call("rbind", meta.list)

pic.list <- list()
for (sample.i in sample.list) {
  pic.list[[sample.i]] <- ggplot(meta.df[meta.df$sample == sample.i, ]) +
    geom_smooth(
      aes(x = loc, y = num, group = type, color = type),
      method = "glm", se = FALSE, formula = y ~ poly(x, 3),
      method.args = list(family = quasibinomial)
    ) + theme_classic() +
    scale_color_manual(values = cell.type.colors) +
    scale_x_continuous(
      breaks = c(0, 45, 90, 135),
      labels = c("0°", "45°", "90°", "135°")
    ) +
    scale_y_continuous(
      breaks = seq(0, 10) / 10,
      labels = seq(0, 10) * 10
    ) +
    coord_cartesian(ylim = c(0, 0.4)) +
    geom_vline(xintercept = c(45, 90), linetype = 3) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 0.1),
      axis.text = element_text(color = "black")
    )
}

pic.list[["PCW9"]] | pic.list[["PCW12_early"]] | pic.list[["PCW12_late"]] |
  pic.list[["PCW14"]] | pic.list[["PCW16"]] | pic.list[["PCW17"]]

ggsave(
  filename = "outputs/Visualization/fig-S5c.pdf",
  width = 20, height = 4
)
