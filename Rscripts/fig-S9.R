# Packages & settings ####
library(Seurat)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Deconvolution/decon_results.RData")
load("outputs/WGCNA/sft.RData")
load("outputs/WGCNA/net.RData")
load("outputs/DEG/me_go.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)
retina.st.domain$module <- paste0(
  "ME",
  apply(
    net$MEs[Cells(retina.st.domain), paste0("ME", 1:4)],
    1, which.max
  )
)

# Scale independence & mean connectivity (Fig. S9a-S9b) ####
## Scale independence ####
plot.data <- data.frame(
  x = sft$fitIndices$Power,
  y = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
  label = seq(1, nrow(sft$fitIndices))
)

pic.a <- ggplot(plot.data) +
  geom_text(aes(x = x, y = y, label = label), color = "red") +
  geom_hline(yintercept = 0.85, linetype = 3, color = "red") +
  xlab("Soft threshold (power)") +
  ylab("Scale-free topology model fit, signed R^2") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black")
  )

## Mean connectivity ####
plot.data <- data.frame(
  x = sft$fitIndices$Power,
  y = sft$fitIndices$mean.k.,
  label = seq(1, nrow(sft$fitIndices))
)

pic.b <- ggplot(plot.data) +
  geom_text(aes(x = x, y = y, label = label), color = "red") +
  geom_hline(yintercept = 0, linetype = 3, color = "red") +
  xlab("Soft threshold (power)") +
  ylab("Mean connectivity") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black")
  )

## Subplots concatenation ####
pic.a + pic.b

ggsave(
  "outputs/Visualization/fig-S9ab.pdf",
  width = 8, height = 4
)

# Weight adjacency heat map (Fig. S9c) ####
me.data <- net$MEs
me.data <- me.data[, -ncol(me.data)]

pdf(file = "outputs/Visualization/fig-S9c.pdf", width = 5, height = 4)
plotEigengeneNetworks(
  me.data, setLabels = "Weight adjacency heatmap",
  plotDendrograms = FALSE, xLabelsAngle = 0, marHeatmap = c(3, 3, 3, 3)
)
dev.off()

# Spatial distribution for WGCNA (Fig. S9d) ####
img.list <- img.list[1:5]
coo.list <- coo.list[1:5]

max.y <- dim(img.list[[1]])[1]
start.x <- 0
for (i in 2:length(img.list)) {
  start.x <- c(start.x, start.x[i - 1] + dim(img.list[[i - 1]])[2] + 100)
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

coo$module <- retina.st.domain@meta.data[rownames(coo), "module"]
pic +
  geom_point(
    data = coo,
    aes(x = coord.x, y = coord.y, fill = module),
    size = 1, shape = 21, color = "white", stroke = 0.1
  ) +
  scale_fill_manual(values = module.colors) +
  coord_fixed(
    xlim = c(0, start.x[5] + dim(img.list[[5]])[2]),
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
  "outputs/Visualization/fig-S9d.pdf",
  width = 12, height = 4
)

# Gene module domain enrichment (Fig. S9e) ####
shuffle.idx <- sample(1:ncol(retina.st.domain))
spot.idx <- order(retina.st.domain$domain[shuffle.idx])
anno_col <- data.frame(
  Domain = retina.st.domain$domain[shuffle.idx][spot.idx],
  row.names = Cells(retina.st.domain)[shuffle.idx][spot.idx]
)
anno_colors <- list(Domain = domain.colors)

bk <- c(seq(-1, -0.01, by = 0.01), seq(0, 1, by = 0.01))
pheatmap(
  t(net$MEs)[paste0("ME", 1:4), shuffle.idx[spot.idx]], scale = "row",
  cluster_rows = FALSE, cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE, show_colnames = FALSE,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  color = c(
    colorRampPalette(colors = c("#5D94A4", "white"))(length(bk) / 2),
    colorRampPalette(colors = c("white", "#DA3B46"))(length(bk) / 2)
  ), breaks = bk,
  filename = "outputs/Visualization/fig-S9e.pdf",
  width = 6, height = 4.5
)

# Gene module GO (Fig. S9f) ####
dotplot(
  me.go,
  label_format = 100,
  font.size = 9,
  showCategory = 7
)

ggsave(
  "outputs/Visualization/fig-S9f.pdf",
  width = 7, height = 6
)
