# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggforce)
library(ggnewscale)
library(VennDiagram)
library(pheatmap)
library(grid)
library(viridis)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/SC/retina_sc_filter.RData")
load("outputs/SC/cb_sc_filter.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/DEG/consensus_deg_d4_d8.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")

# D4 & D8 DEG intersection (Fig. S8a) ####
domain.deg <- domain.deg[domain.deg$avg_log2FC > 0,]
venn.list <- list(
  D4 = unique(domain.deg$gene[domain.deg$cluster == "D4"]),
  D8 = unique(domain.deg$gene[domain.deg$cluster == "D8"])
)
venn.diagram(
  venn.list, disable.logging = TRUE,
  imagetype = "svg", cex = 2, cat.cex = 2,
  fill = c("#66C2A5", "#FC8D62"),
  alpha = c(0.6, 0.6), lwd = 3,
  filename = "outputs/Visualization/fig-S8a.svg",
  width = 6, height = 6
)

# Consensus DEG spatial expression (Fig. S8b) ####
max.y <- dim(img.list[[1]])[1]
start.x <- 0
for (i in 2:length(img.list)) {
  start.x <- c(start.x, start.x[i - 1] + dim(img.list[[i - 1]])[2] + 150)
  coo.list[[i]]$coord.x <- coo.list[[i]]$coord.x + start.x[i]
  if (max.y < dim(img.list[[i]])[1]) {
    max.y <- dim(img.list[[i]])[1]
  }
}
start.y <- NULL
for (i in seq_along(coo.list)) {
  start.y <- c(start.y, ceiling((max.y - dim(img.list[[i]])[1]) / 2))
  coo.list[[i]]$coord.y <- coo.list[[i]]$coord.y + start.y[i]
}
coo <- do.call("rbind", coo.list)
rownames(coo) <- sub("^[^.]*\\.", "", rownames(coo))
coo$exp <- colMeans(retina.st.domain@assays$Spatial@data[co.deg.d4.d8, rownames(coo)])
coo$exp[coo$exp > quantile(coo$exp, 0.95)] <- quantile(coo$exp, 0.95)

pic <- ggplot()
for (i in seq_along(img.list)) {
  img.i <- img.list[[i]]
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
  geom_point(
    data = coo,
    aes(x = coord.x, y = coord.y, fill = exp),
    size = 0.65, shape = 21, stroke = 0.05, color = "white"
  ) +
  scale_fill_viridis(option = "B") +
  coord_fixed(
    xlim = c(0, start.x[6] + dim(img.list[[6]])[2]),
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
  "outputs/Visualization/fig-S8b.pdf",
  width = 10, height = 5
)

# Consensus DEG retina single-cell UMAP (Fig. S8c-S8d) ####
common.gene <- intersect(co.deg.d4.d8, rownames(retina.sc.filter))
retina.sc.filter$avg.exp <- colMeans(
  retina.sc.filter@assays$RNA@data[common.gene,]
)

pic.c <- DimPlot(
  retina.sc.filter, pt.size = 2,
  group.by = "cell.type", raster = TRUE,
  cols = cell.type.colors, label = TRUE
) +
  theme_bw() +
  coord_fixed(ratio = 0.85) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
pic.d <- FeaturePlot(
  retina.sc.filter, pt.size = 2,
  features = "avg.exp", raster = TRUE
) +
  scale_color_viridis(option = "B") +
  theme_bw() +
  coord_fixed(ratio = 0.85) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

pic.c | pic.d
ggsave(
  filename = "outputs/Visualization/fig-S8cd.pdf",
  width = 10, height = 3.5
)

# Consensus DEG retina single-cell expression heat map (Fig. S8e) ####
co.deg.d4.d8.ava <- intersect(co.deg.d4.d8, rownames(retina.sc.filter))
exp.mat <- retina.sc.filter@assays$RNA@data[co.deg.d4.d8.ava,]
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
  scale_fill_viridis(option = "B") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Avg. expr.")

ggsave(
  filename = "outputs/Visualization/fig-S8e.pdf",
  width = 4, height = 4
)

# Consensus DEG ciliary body single-cell UMAP (Fig. S8f-S8g) ####
common.gene <- intersect(co.deg.d4.d8, rownames(cb.sc.filter))
cb.sc.filter$avg.exp <- colMeans(cb.sc.filter@assays$RNA@data[common.gene,])

pic.f <- DimPlot(
  cb.sc.filter, pt.size = 2,
  group.by = "cell.type", raster = TRUE, label = TRUE
) +
  theme_bw() +
  coord_fixed(ratio = 0.85) +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )
pic.g <- FeaturePlot(
  cb.sc.filter, pt.size = 2,
  features = "avg.exp", raster = TRUE
) +
  scale_color_viridis(option = "B") +
  theme_bw() +
  coord_fixed(ratio = 0.85) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_blank()
  )

pic.f | pic.g
ggsave(
  filename = "outputs/Visualization/fig-S8fg.pdf",
  width = 10, height = 3.5
)

# Consensus DEG ciliary body single-cell expression heat map (Fig. S8h) ####
co.deg.d4.d8.ava <- intersect(co.deg.d4.d8, rownames(cb.sc.filter))
exp.mat <- cb.sc.filter@assays$RNA@data[co.deg.d4.d8.ava,]
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
  scale_fill_viridis(option = "B") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Avg. expr.")

ggsave(
  filename = "outputs/Visualization/fig-S8h.pdf",
  width = 4, height = 4
)
