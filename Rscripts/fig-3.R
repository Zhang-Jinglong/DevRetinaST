# Packages & settings ####
library(Seurat)
library(ggplot2)
library(pheatmap)
library(ggforce)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Deconvolution/decon_results.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")
load("outputs/ImageProcessing/coordinate.RData")
load("outputs/WGCNA/net.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$layer <- radial.anno[spot.list, "layer"]
retina.st.domain$quantile <- radial.anno[spot.list, "quantile"]
retina.st.domain$tangential <- degree.df[spot.list, "degree"]
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# UMAP for clusters (Fig. 3a) ####
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2],
  domain = retina.st.domain$domain
)
ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = domain
  )) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(values = domain.colors) +
  coord_fixed(ratio = 2)

ggsave(
  "outputs/Visualization/fig-3a.pdf",
  width = 5, height = 5
)

# Cluster enrich cell type fraction (Fig. 3b) ####
domain.type.list <- list()
for (domain.i in unique(retina.st.domain$domain)) {
  idx.i <- Cells(retina.st.domain)[retina.st.domain$domain == domain.i]
  domain.type.list[[domain.i]] <- rowMeans(pred.cell2loction[, idx.i])
}
domain.type <- do.call("rbind", domain.type.list)

domain.type <- domain.type[domain.list, ]
domain.type <- domain.type[, cell.type.list[c(6, 3, 5, 1, 2, 4, 10, 8, 9, 7)]]

bk <- c(seq(-3, -0.01, by = 0.01), seq(0, 3, by = 0.01))
pheatmap(
  t(domain.type), scale = "row",
  cluster_cols = FALSE, cluster_rows = FALSE,
  breaks = bk, angle_col = "45", fontsize = 20,
  color = c(
    colorRampPalette(colors = c("#5D94A4", "white"))(length(bk) / 2),
    colorRampPalette(colors = c("white", "#DA3B46"))(length(bk) / 2)
  ),
  border_color = "#000000",
  filename = "outputs/Visualization/fig-3b.pdf",
  width = 9.5, height = 6
)
dev.new()

# Radial/Tangential cluster difference (Fig. 3c) ####
meta <- retina.st.domain@meta.data

## Radial difference ####
meta$quantile[meta$layer == "GCL"] <- meta$quantile[meta$layer == "GCL"] + 1
meta$quantile <- 2 - meta$quantile

order.i <- aggregate(
  meta$quantile, by = list(group = meta$domain), median
)
order.i <- order.i$group[order(order.i$x, decreasing = TRUE)]
meta$domain <- factor(meta$domain, levels = order.i)
meta$name <- "Radial layer"

pic.c <- ggplot(meta) +
  geom_boxplot(
    aes(y = quantile, x = domain, fill = domain),
    outlier.alpha = 0, size = 0.1
  ) +
  scale_fill_manual(values = domain.colors) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0),
    labels = c("Border", radial.list[6:1])
  ) +
  theme_bw() + ylab("Layer") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  ) + facet_grid(. ~ name)

## Tangential difference ####
meta$tangential.abs <- abs(meta$tangential)

order.i <- aggregate(
  meta$tangential.abs, by = list(group = meta$domain), median
)
order.i <- order.i$group[order(order.i$x, decreasing = TRUE)]
meta$domain <- factor(meta$domain, levels = order.i)
meta$name <- "Absolute value of tangential angle"

pic.d <- ggplot(meta) +
  geom_boxplot(
    aes(y = tangential.abs, x = domain, fill = domain),
    outlier.alpha = 0, size = 0.1
  ) +
  scale_fill_manual(values = domain.colors) +
  scale_y_continuous(
    breaks = c(0, 45, 90, 135),
    labels = c("0°", "45°", "90°", "135°")
  ) +
  theme_bw() + ylab("Angle") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  ) + facet_grid(. ~ name)

## Subplots concatenation ####
pic.c / pic.d

ggsave(
  "outputs/Visualization/fig-3c.pdf",
  width = 4, height = 4
)

# Spatial distribution for clusters (Fig. 3d) ####
plot.data <- data.frame(
  value = retina.st.domain$domain,
  x = (
    retina.st.domain$coord.x +
      coord.add.1[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.1[retina.st.domain$sample.name, "y"]
  )
)

pic.d <- ggplot()
for (sample.i in sample.list) {
  pic.d <- create.empty.bk(
    pic.d, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.d +
  geom_point(
    data = plot.data, size = 0.4,
    mapping = aes(x = x, y = y, color = value)
  ) +
  scale_color_manual(values = domain.colors) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  coord_fixed()

ggsave(
  "outputs/Visualization/fig-3d.pdf",
  width = 10, height = 4
)

# Spatial WGCNA (Fig. 3e) ####
retina.st.domain$module <- paste0(
  "ME", 
  apply(
    net$MEs[Cells(retina.st.domain), paste0("ME", 1:4)],
    1, which.max
  )
)

plot.data <- data.frame(
  value = retina.st.domain$module,
  x = (
    retina.st.domain$coord.x +
      coord.add.1[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.1[retina.st.domain$sample.name, "y"]
  )
)

pic.e <- ggplot()
for (sample.i in sample.list) {
  pic.e <- create.empty.bk(
    pic.e, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.e +
  geom_point(
    data = plot.data, size = 0.4,
    mapping = aes(x = x, y = y, color = value)
  ) +
  scale_color_manual(values = module.colors) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  ) +
  coord_fixed()

ggsave(
  "outputs/Visualization/fig-3e.pdf",
  width = 10, height = 4
)

# Box plot of cluster module (Fig. 3f) ####
plot.list <- list()
for (me.i in paste0("ME", 1:4)) {
  df <- data.frame(
    value = net$MEs[Cells(retina.st.domain), me.i],
    group = retina.st.domain$domain,
    me = me.i
  )
  order.i <- aggregate(df$value, by = list(group = df$group), median)
  order.i <- order.i$group[order(order.i$x, decreasing = TRUE)]
  df$group <- factor(df$group, levels = order.i)
  plot.list[[me.i]] <- ggplot(
    df, aes(x = group, y = value, color = group)
  ) +
    geom_boxplot(outlier.alpha = 0, size = 0.1) +
    scale_y_continuous(breaks = c(-5:5 * 0.04)) +
    scale_color_manual(values = domain.colors) +
    ylab("Eigengene score") +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black")
    ) + facet_grid(. ~ me)
}

plot.list[["ME1"]] | plot.list[["ME2"]] |
  plot.list[["ME3"]] | plot.list[["ME4"]]

ggsave(
  filename = "outputs/Visualization/fig-3f.pdf",
  width = 11, height = 2
)
