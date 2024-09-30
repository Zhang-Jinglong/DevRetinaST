# Packages & settings ####
library(Seurat)
library(ggplot2)
library(grid)
library(patchwork)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")

# Example markers on H&E (Fig. 1c) ####
lim.1 <- c(28, 78, 93, 143)
lim.2 <- c(70, 120, 240, 290)
lim.3 <- c(8, 58, 210, 260)

img.sub.list <- img.list[c(1, 4, 6)]
max.y <- dim(img.sub.list[[1]])[1]
start.x <- 0
for (i in 2:length(img.sub.list)) {
  start.x <- c(start.x, start.x[i - 1] +
    dim(img.sub.list[[i - 1]])[2] +
    50)
  if (max.y < dim(img.sub.list[[i]])[1]) {
    max.y <- dim(img.sub.list[[i]])[1]
  }
}
start.y <- NULL
for (i in seq_along(img.sub.list)) {
  start.y <- c(start.y, ceiling((max.y - dim(img.sub.list[[i]])[1]) / 2))
}

pic <- ggplot()
for (i in seq_along(img.sub.list)) {
  img.i <- img.sub.list[[i]]
  pic <- pic +
    annotation_custom(
      rasterGrob(
        img.i,
        width = unit(1, "npc"),
        height = unit(1, "npc")
      ),
      xmin = start.x[i], xmax = start.x[i] + dim(img.i)[2],
      ymin = start.y[i], ymax = start.y[i] + dim(img.i)[1]
    )
}
pic +
  annotate(
    "rect",
    xmin = lim.1[1] + start.x[1], xmax = lim.1[2] + start.x[1],
    ymin = lim.1[3] + start.y[1], ymax = lim.1[4] + start.y[1],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "rect",
    xmin = lim.2[1] + start.x[2], xmax = lim.2[2] + start.x[2],
    ymin = lim.2[3] + start.y[2], ymax = lim.2[4] + start.y[2],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  annotate(
    "rect",
    xmin = lim.3[1] + start.x[3], xmax = lim.3[2] + start.x[3],
    ymin = lim.3[3] + start.y[3], ymax = lim.3[4] + start.y[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  coord_fixed(
    xlim = c(0, start.x[3] + dim(img.sub.list[[3]])[2]),
    ylim = c(0, max.y), expand = FALSE
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.spacing = element_blank()
  )

ggsave(
  "outputs/Visualization/fig-1c_background.pdf",
  width = 18, height = 10
)

zoom.plot <- function(sample.name, gene.name, lim.i, pt.size, max.value) {
  coo <- coo.list[[sample.name]]
  img <- img.list[[sample.name]]

  coo$exp <- retina.st.domain@assays$Spatial@data[gene.name, rownames(coo)]
  coo$exp[coo$exp > quantile(coo$exp, 0.9)] <- quantile(coo$exp, 0.9)
  lim.0 <- c(0, dim(img)[2], 0, dim(img)[1])

  pic.i <- ggplot(
    coo, aes(x = coord.x - lim.i[1], y = coord.y - lim.i[3], fill = exp)
  ) +
    annotation_custom(
      rasterGrob(
        img[(lim.0[4] - lim.i[4] + 1):(lim.0[4] - lim.i[3]), (lim.i[1] + 1):lim.i[2],],
        width = unit(1, "npc"),
        height = unit(1, "npc")
      ),
      xmin = 0, xmax = lim.i[2] - lim.i[1],
      ymin = 0, ymax = lim.i[4] - lim.i[3]
    ) +
    geom_point(
      size = pt.size, shape = 21, color = "white", stroke = pt.size * 0.1
    ) +
    coord_fixed(
      xlim = c(0, lim.i[2] - lim.i[1]),
      ylim = c(0, lim.i[4] - lim.i[3]),
      expand = FALSE
    ) +
    scale_fill_viridis(option = "B", limits = c(0, max.value)) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.spacing = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "pt")
    )
  pic.i
}

pic.d1 <- zoom.plot("PCW9", "SOX2", lim.1, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.1[2] - lim.1[1], ymin = 0, ymax = lim.1[4] - lim.1[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d2 <- zoom.plot("PCW9", "SOX9", lim.1, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.1[2] - lim.1[1], ymin = 0, ymax = lim.1[4] - lim.1[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d3 <- zoom.plot("PCW9", "TUBB3", lim.1, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.1[2] - lim.1[1], ymin = 0, ymax = lim.1[4] - lim.1[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d4 <- zoom.plot("PCW9", "SNCG", lim.1, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.1[2] - lim.1[1], ymin = 0, ymax = lim.1[4] - lim.1[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")

pic.d5 <- zoom.plot("PCW14", "SOX2", lim.2, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.2[2] - lim.2[1], ymin = 0, ymax = lim.2[4] - lim.2[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d6 <- zoom.plot("PCW14", "SOX9", lim.2, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.2[2] - lim.2[1], ymin = 0, ymax = lim.2[4] - lim.2[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d7 <- zoom.plot("PCW14", "TUBB3", lim.2, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.2[2] - lim.2[1], ymin = 0, ymax = lim.2[4] - lim.2[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d8 <- zoom.plot("PCW14", "SNCG", lim.2, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.2[2] - lim.2[1], ymin = 0, ymax = lim.2[4] - lim.2[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")

pic.d9 <- zoom.plot("PCW17", "SOX2", lim.3, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.3[2] - lim.3[1], ymin = 0, ymax = lim.3[4] - lim.3[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d10 <- zoom.plot("PCW17", "SOX9", lim.3, 7, 2) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.3[2] - lim.3[1], ymin = 0, ymax = lim.3[4] - lim.3[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d11 <- zoom.plot("PCW17", "TUBB3", lim.3, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.3[2] - lim.3[1], ymin = 0, ymax = lim.3[4] - lim.3[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")
pic.d12 <- zoom.plot("PCW17", "SNCG", lim.3, 7, 4.5) +
  annotate(
    "rect",
    xmin = 0, xmax = lim.3[2] - lim.3[1], ymin = 0, ymax = lim.3[4] - lim.3[3],
    alpha = 0, color = "black", linetype = "dashed", linewidth = 1
  ) +
  theme(legend.position = "none")

(pic.d1 +
  pic.d2 +
  pic.d3 +
  pic.d4 +
  plot_layout(nrow = 2)) |
  (pic.d5 +
    pic.d6 +
    pic.d7 +
    pic.d8 +
    plot_layout(nrow = 2)) |
  (pic.d9 +
    pic.d10 +
    pic.d11 +
    pic.d12 +
    plot_layout(nrow = 2))

ggsave(
  "outputs/Visualization/fig-1c_zoom.pdf",
  width = 15, height = 5
)

zoom.plot("PCW17", "SOX2", lim.3, 7, 2) +
  zoom.plot("PCW17", "SOX9", lim.3, 7, 2) +
  zoom.plot("PCW17", "TUBB3", lim.3, 7, 4.5) +
  zoom.plot("PCW17", "SNCG", lim.3, 7, 4.5)

ggsave(
  "outputs/Visualization/fig-1c_legend.pdf",
  width = 5, height = 5
)

# UMAP for all spots (Fig. 1d) ####
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2],
  radial = radial.anno[Cells(retina.st.domain), "sublayer"],
  tangential = abs(degree.df[Cells(retina.st.domain), "degree"])
)
umap.data$radial <- factor(umap.data$radial, levels = radial.list)

ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = radial,
    size = tangential
  ), alpha = 0.8, shape = 16) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_manual(values = radial.colors) +
  scale_size_continuous(breaks = c(0, 45, 90, 135), range = c(5, 0.5)) +
  coord_fixed(ratio = 2)

ggsave(
  "outputs/Visualization/fig-1d.pdf",
  width = 7, height = 7
)
