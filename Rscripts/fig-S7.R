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

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Scale independence & mean connectivity (Fig. S7a-S7b) ####
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
  "outputs/Visualization/fig-S7ab.pdf",
  width = 8, height = 4
)

# Weight adjacency heat map (Fig. S7c) ####
me.data <- net$MEs
me.data <- me.data[, -ncol(me.data)]

pdf(file = "outputs/Visualization/fig-S7c.pdf", width = 5, height = 4)
plotEigengeneNetworks(
  me.data, setLabels = "Weight adjacency heatmap",
  plotDendrograms = FALSE, xLabelsAngle = 0, marHeatmap = c(3, 3, 3, 3)
)
dev.off()

# Gene module domain enrichment (Fig. S7d) ####
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
  filename = "outputs/Visualization/fig-S7d.pdf",
  width = 6, height = 4.5
)

# Gene module GO (Fig. S7e) ####
dotplot(
  me.go,
  label_format = 100,
  font.size = 9,
  showCategory = 7
)

ggsave(
  "outputs/Visualization/fig-S7e.pdf",
  width = 7, height = 6
)
