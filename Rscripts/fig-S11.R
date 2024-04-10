# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggplot2)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
## Spatial data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

## LR list ####
lr.table <- read.csv(
  "outputs/CCC/COMMOT/retina_db_filter.csv", row.names = 1
)
rownames(lr.table) <- paste(lr.table$X0, lr.table$X1, sep = "-")

lr.summary.list <- list()
for (sample.i in sample.list) {
  lr.summary.i <- read.csv(
    paste0("outputs/CCC/COMMOT/commot_", sample.i, "_lr.csv"),
    header = FALSE
  )
  colnames(lr.summary.i) <- c("LR", "value")
  lr.summary.i$sample <- sample.i
  lr.summary.i$value <- as.numeric(scale(lr.summary.i$value))
  
  lr.summary.list[[sample.i]] <- lr.summary.i
}
lr.summary <- do.call("rbind", lr.summary.list)
lr.summary.df <- aggregate(
  lr.summary$value, by = list(LR = lr.summary$LR), mean
)
lr.summary.df$type <- lr.table[lr.summary.df$LR, "X3"]
lr.target <- lr.summary.df$LR[lr.summary.df$x > 2]
rm(lr.summary, lr.summary.df, lr.summary.list, lr.table, lr.summary.i)

## CCC results ####
cpdb.means <- read.csv(
  "outputs/CCC/CellPhoneDB/cpdb_filter_means.csv",
  row.names = 1
)
rownames(cpdb.means) <- paste(cpdb.means$gene_a, cpdb.means$gene_b, sep = "-")
cpdb.means <- cpdb.means[rownames(cpdb.means) %in% lr.target, ]
cpdb.pvals <- read.csv(
  "outputs/CCC/CellPhoneDB/cpdb_filter_pvalues.csv",
  row.names = 1
)
rownames(cpdb.pvals) <- paste(cpdb.pvals$gene_a, cpdb.pvals$gene_b, sep = "-")
cpdb.pvals <- cpdb.pvals[rownames(cpdb.pvals) %in% lr.target, ]

cpdb.pvals.sub <- cpdb.pvals[, 3:102]
cpdb.ava.idx <- colnames(cpdb.pvals.sub)[colSums(cpdb.pvals.sub < 0.05) > 0]
cpdb.means <- cpdb.means[, cpdb.ava.idx]
cpdb.pvals <- cpdb.pvals[, cpdb.ava.idx]
rm(cpdb.pvals.sub, cpdb.ava.idx)

# CCC dot plot (Fig. S11b) ####
cpdb.means.df <- reshape2::melt(as.matrix(cpdb.means))
cpdb.pvals.df <- reshape2::melt(as.matrix(cpdb.pvals))
cpdb.means.df$pvals <- - log10(cpdb.pvals.df$value + 1e-3)

cpdb.max <- quantile(cpdb.means.df$value, 0.95)
cpdb.min <- quantile(cpdb.means.df$value, 0.05)
cpdb.means.df$value[cpdb.means.df$value > cpdb.max] <- cpdb.max
cpdb.means.df$value[cpdb.means.df$value < cpdb.min] <- cpdb.min
cpdb.means.df$value <- scale(cpdb.means.df$value)

ggplot(cpdb.means.df) +
  geom_point(
    aes(x = Var2, y = Var1, color = value, size = pvals)
  ) +
  scale_size_continuous(range = c(0, 4)) +
  scale_color_gradient2(
    low = "#5D94A4", mid = "#CCCCCC", high = "#DA3B46"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(color = "Normalized CCI means", size = "P value (-log10)")

ggsave(
  "outputs/Visualization/fig-S11b.pdf",
  width = 12, height = 3.6
)

# CNTN1 & NOTCH1 co-localization (Fig. S11c) ####
plot.data <- data.frame(
  l.value = retina.st.domain@assays$Spatial@data["CNTN1", ],
  r.value = retina.st.domain@assays$Spatial@data["NOTCH1", ],
  x = (
    retina.st.domain$coord.x +
      coord.add.1[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.1[retina.st.domain$sample.name, "y"]
  )
)
plot.data$l.value <- scale(plot.data$l.value)
plot.data$l.value[plot.data$l.value > 2.5] <- 2.5
plot.data$l.value[plot.data$l.value < -2.5] <- -2.5
plot.data$r.value <- scale(plot.data$r.value)
plot.data$r.value[plot.data$r.value > 2.5] <- 2.5
plot.data$r.value[plot.data$r.value < -2.5] <- -2.5

pic.c <- ggplot()
for (sample.i in sample.list) {
  pic.c <- create.empty.bk(
    pic.c, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.c +
  geom_point(
    data = plot.data, size = 0.5,
    mapping = aes(x = x, y = y, color = l.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#EE692E",
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_color() +
  geom_point(
    data = plot.data, size = 0.1,
    mapping = aes(x = x, y = y, color = r.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#9A1DEE",
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_color() +
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
  "outputs/Visualization/fig-S11c.pdf",
  width = 10, height = 4
)

# PTN & PTPRS co-localization (Fig. S11d) ####
plot.data <- data.frame(
  l.value = retina.st.domain@assays$Spatial@data["PTN", ],
  r.value = retina.st.domain@assays$Spatial@data["PTPRS", ],
  x = (
    retina.st.domain$coord.x +
      coord.add.1[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.1[retina.st.domain$sample.name, "y"]
  )
)
plot.data$l.value <- scale(plot.data$l.value)
plot.data$l.value[plot.data$l.value > 2.5] <- 2.5
plot.data$l.value[plot.data$l.value < -2.5] <- -2.5
plot.data$r.value <- scale(plot.data$r.value)
plot.data$r.value[plot.data$r.value > 2.5] <- 2.5
plot.data$r.value[plot.data$r.value < -2.5] <- -2.5

pic.d <- ggplot()
for (sample.i in sample.list) {
  pic.d <- create.empty.bk(
    pic.d, coord.add.1[sample.i, "x"], coord.add.1[sample.i, "y"]
  )
}

pic.d +
  geom_point(
    data = plot.data, size = 0.5,
    mapping = aes(x = x, y = y, color = l.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#EE692E",
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_color() +
  geom_point(
    data = plot.data, size = 0.1,
    mapping = aes(x = x, y = y, color = r.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#9A1DEE",
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_color() +
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
  "outputs/Visualization/fig-S11d.pdf",
  width = 10, height = 4
)
