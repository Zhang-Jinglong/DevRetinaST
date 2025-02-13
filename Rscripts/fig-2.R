# Packages & settings ####
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(scales)
library(ggforce)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")
load("outputs/QualityControl/retina_sc_filter.RData")
load("outputs/Deconvolution/decon_results.RData")
load("outputs/ImageProcessing/radial_anno.RData")
load("outputs/ImageProcessing/tangential_anno.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.filter)
retina.st.filter$radial <- radial.anno[spot.list, "sublayer"]
retina.st.filter$tangential <- degree.df[spot.list, "degree"]
retina.st.filter$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.filter$coord.y <- spatial.coord.new[spot.list, "coord.y"]

# Gene correlation between single-cell & Visium (Fig. 2a) ####
retina.sc.scale <- retina.sc.qc %>%
  NormalizeData(verbose = FALSE) %>%
  ScaleData(features = rownames(retina.sc.qc))
retina.sc.split <- SplitObject(retina.sc.scale, split.by = "sample.name")

retina.st.scale <- retina.st.filter %>%
  NormalizeData(verbose = FALSE) %>%
  ScaleData(features = rownames(retina.st.filter))
retina.st.split <- SplitObject(retina.st.scale, split.by = "sample.name")

corr.mat <- matrix(
  0, nrow = length(retina.sc.split), ncol = length(retina.st.split)
)
rownames(corr.mat) <- names(retina.sc.split)
colnames(corr.mat) <- names(retina.st.split)
for (st.i in names(retina.st.split)) {
  for (sc.i in names(retina.sc.split)) {
    retina.st.i <- retina.st.split[[st.i]]@assays$Spatial@scale.data
    retina.sc.i <- retina.sc.split[[sc.i]]@assays$RNA@scale.data
    retina.st.mat <- as.matrix(retina.st.i)
    retina.sc.mat <- as.matrix(retina.sc.i)
    st.gene.exp.median <- rowMeans(retina.st.mat)
    sc.gene.exp.median <- rowMeans(retina.sc.mat)
    gene.list <- intersect(
      names(st.gene.exp.median), names(sc.gene.exp.median)
    )
    st.gene.exp.median <- st.gene.exp.median[gene.list]
    sc.gene.exp.median <- sc.gene.exp.median[gene.list]
    corr.mat[sc.i, st.i] <- cor(
      sc.gene.exp.median, st.gene.exp.median, method = "spearman"
    )
  }
}
corr.mat <- corr.mat[c(8, 1:7), 6:1]

bk <- c(seq(-3, -0.01, by = 0.01), seq(0, 3, by = 0.01))
pheatmap(
  t(corr.mat), scale = "none",
  cluster_cols = FALSE, cluster_rows = FALSE,
  angle_col = "45", fontsize = 20,
  color = c(
    colorRampPalette(colors = c("#5D94A4", "white"))(length(bk) / 2),
    colorRampPalette(colors = c("white", "#DA3B46"))(length(bk) / 2)
  ),
  border_color = "#000000",
  filename = "outputs/Visualization/fig-2a.pdf",
  width = 9.5, height = 6
)

# Spatial-temporal cell-type distribution (Fig. 2b) ####
coord.add.t <- data.frame(
  x = rep(0, 6),
  y = seq(0, 5, 1) * (-15),
  row.names = sort(cell.type.list[1:6])
)

pic.c <- ggplot()
for (type.i in sort(cell.type.list[1:6])) {
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
    pic.c <- create.empty.bk(
      pic.c,
      coord.add.t[type.i, "x"] + coord.add.1[sample.i, "x"],
      coord.add.t[type.i, "y"] + coord.add.1[sample.i, "y"]
    )
  }

  pic.c <- pic.c +
    geom_point(
      data = plot.data.i, size = 0.4,
      mapping = aes(x = x, y = y, color = value)
    ) +
    scale_color_gradient(low = "#FFFFFF", high = cell.type.colors[type.i]) +
    new_scale_color()
}
pic.c +
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
  "outputs/Visualization/fig-2b.pdf",
  width = 10, height = 8
)

# Pseudo cell type proportion change curve (Fig. 2c) ####
type.curve.change <- list()
for (sample.i in sample.list) {
  type.curve <- list(
    sample = c(), cell.type = c(), sum.frac = c()
  )
  for (type.i in cell.type.list) {
    sample.idx <- retina.st.filter$sample.name == sample.i
    sum.frac.i <- sum(
      correct.cell2loction[type.i, Cells(retina.st.filter)[sample.idx]]
    )
    type.curve$sample <- c(type.curve$sample, sample.i)
    type.curve$cell.type <- c(type.curve$cell.type, type.i)
    type.curve$sum.frac <- c(type.curve$sum.frac, sum.frac.i)
  }

  type.curve <- as.data.frame(type.curve)
  type.curve$inferred.frac <- (
    type.curve$sum.frac / sum(type.curve$sum.frac)
  )
  type.curve.change[[sample.i]] <- type.curve
}

type.curve.change <- do.call(rbind, type.curve.change)
type.curve.change$sample <- factor(
  type.curve.change$sample, levels = sample.list
)
type.curve.change$inferred.frac <- type.curve.change$inferred.frac * 100
type.curve.change$group <- "Precursor cell types"
diff.idx <- type.curve.change$cell.type %in% cell.type.list[1:6]
type.curve.change$group[diff.idx] <- "Differentiated cell types"

pic.c1 <- ggplot(
  type.curve.change[type.curve.change$group == "Differentiated cell types",],
  aes(
    x = sample, y = inferred.frac,
    color = cell.type, group = cell.type
  )
) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 18) +
  scale_y_continuous(limits = c(0, 50)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "black")
  ) +
  scale_color_manual(values = cell.type.colors[1:6]) +
  facet_grid(. ~ group)

pic.c2 <- ggplot(
  type.curve.change[type.curve.change$group != "Differentiated cell types",],
  aes(
    x = sample, y = inferred.frac,
    color = cell.type, group = cell.type
  )
) +
  geom_line(linewidth = 0.3) +
  geom_point(shape = 18) +
  scale_y_continuous(limits = c(0, 50)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(color = "black")
  ) +
  scale_color_manual(values = cell.type.colors[7:10]) +
  facet_grid(. ~ group)

pic.c1 | pic.c2
ggsave(
  "outputs/Visualization/fig-2c.pdf",
  width = 7, height = 4
)
