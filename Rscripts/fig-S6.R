# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(VennDiagram)
library(pheatmap)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/SC/retina_sc_filter.RData")
load("outputs/SC/cb_sc_filter.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/DEG/consensus_deg_d4_d8.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

# D4 & D8 DEG intersection (Fig. S6a) ####
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0, ]
venn.list <- list(
  D4 = unique(domain.deg$gene[domain.deg$cluster == "D4"]),
  D8 = unique(domain.deg$gene[domain.deg$cluster == "D8"])
)
venn.diagram(
  venn.list, disable.logging = TRUE,
  imagetype = "svg", cex = 2, cat.cex = 2,
  fill = c("#66C2A5", "#FC8D62"),
  alpha = c(0.6, 0.6), lwd = 3,
  filename = "outputs/Visualization/fig-S6a.svg",
  width = 6, height = 6
)

# Consensus DEG spatial expression (Fig. S6b)  ####
plot.data <- data.frame(
  value = colMeans(retina.st.domain@assays$Spatial@data[co.deg.d4.d8, ]),
  x = (
    retina.st.domain$coord.x +
      coord.add.2[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.2[retina.st.domain$sample.name, "y"]
  )
)
plot.data$value <- scale(plot.data$value)

pic.c <- ggplot()
for (sample.i in sample.list) {
  pic.c <- create.empty.bk(
    pic.c, coord.add.2[sample.i, "x"], coord.add.2[sample.i, "y"]
  )
}

pic.c +
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
  "outputs/Visualization/fig-S6b.pdf",
  width = 6, height = 5
)

# Consensus DEG retina single-cell UMAP (Fig. S6c-S6d) ####
common.gene <- intersect(co.deg.d4.d8, rownames(retina.sc.filter))
retina.sc.filter$avg.exp <- colMeans(
  retina.sc.filter@assays$RNA@data[common.gene, ]
)
retina.sc.filter$avg.exp <- scale(retina.sc.filter$avg.exp)

pic.c <- DimPlot(
  retina.sc.filter, pt.size = 2,
  group.by = "cell.type", raster = TRUE,
  cols = cell.type.colors, label = TRUE
) + theme_bw() + coord_fixed(ratio = 0.85) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
pic.d <- FeaturePlot(
  retina.sc.filter, pt.size = 2,
  features = "avg.exp", raster = TRUE,
  cols = c("#5D94A4", "#FFFFFF", "#DA3B46")
) + theme_bw() + coord_fixed(ratio = 0.85) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

pic.c | pic.d
ggsave(
  filename = "outputs/Visualization/fig-S6cd.pdf",
  width = 10, height = 3.5
)

# Consensus DEG retina single-cell expression heat map (Fig. S6e) ####
co.deg.d4.d8.ava <- intersect(co.deg.d4.d8, rownames(retina.sc.filter))
exp.mat <- retina.sc.filter@assays$RNA@data[co.deg.d4.d8.ava, ]
exp.mat <- matrix(c(
  rowMeans(exp.mat[, retina.sc.filter$cell.type == "RPC"]),
  rowMeans(exp.mat[, retina.sc.filter$cell.type != "RPC"])
), ncol = 2)

plot.data <- data.frame(
  gene = rep(co.deg.d4.d8.ava, 2),
  group = factor(c(
    rep("RPC", length(co.deg.d4.d8.ava)),
    rep("Others", length(co.deg.d4.d8.ava))
  ), levels = c("RPC", "Others")),
  value = c(exp.mat[, 1], exp.mat[, 2])
)
plot.data$value[plot.data$value > 1] <- 1

ggplot(plot.data) +
  geom_tile(aes(x = group, y = gene, fill = value)) +
  scale_fill_gradient(low = "#BFD3E6", high = "#88419D") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Avg. expr.")

ggsave(
  filename = "outputs/Visualization/fig-S6e.pdf",
  width = 4, height = 4
)

# Consensus DEG ciliary body single-cell UMAP (Fig. S6f-S6g) ####
common.gene <- intersect(co.deg.d4.d8, rownames(cb.sc.filter))
cb.sc.filter$avg.exp <- colMeans(cb.sc.filter@assays$RNA@data[common.gene, ])

pic.f <- DimPlot(
  cb.sc.filter, pt.size = 2,
  group.by = "cell.type", raster = TRUE, label = TRUE
) + theme_bw() + coord_fixed(ratio = 0.85) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
pic.g <- FeaturePlot(
  cb.sc.filter, pt.size = 2,
  features = "avg.exp", raster = TRUE,
  cols = c("#5D94A4", "#FFFFFF", "#DA3B46")
) + theme_bw() + coord_fixed(ratio = 0.85) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

pic.f | pic.g
ggsave(
  filename = "outputs/Visualization/fig-S6fg.pdf",
  width = 10, height = 3.5
)

# Consensus DEG ciliary body single-cell expression heat map (Fig. S6h) ####
co.deg.d4.d8.ava <- intersect(co.deg.d4.d8, rownames(cb.sc.filter))
exp.mat <- cb.sc.filter@assays$RNA@data[co.deg.d4.d8.ava, ]
exp.mat <- matrix(c(
  rowMeans(exp.mat[, cb.sc.filter$cell.type == "NPE"]),
  rowMeans(exp.mat[, cb.sc.filter$cell.type != "NPE"])
), ncol = 2)

plot.data <- data.frame(
  gene = rep(co.deg.d4.d8.ava, 2),
  group = factor(c(
    rep("NPE", length(co.deg.d4.d8.ava)),
    rep("Others", length(co.deg.d4.d8.ava))
  ), levels = c("NPE", "Others")),
  value = c(exp.mat[, 1], exp.mat[, 2])
)
plot.data$value[plot.data$value > 1] <- 1

ggplot(plot.data) +
  geom_tile(aes(x = group, y = gene, fill = value)) +
  scale_fill_gradient(low = "#BFD3E6", high = "#88419D") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Avg. expr.")

ggsave(
  filename = "outputs/Visualization/fig-S6h.pdf",
  width = 4, height = 4
)
