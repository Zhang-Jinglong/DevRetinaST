# Packages & settings ####
library(ggplot2)
library(ggrepel)
library(grid)
library(ggforce)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Domain/diff_res_num_list.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/DEG/domain_deg_go.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Domain resolution selection (Fig. S7a) ####
plot.df <- data.frame(
  group = c(rep("Up-regulated", 11), rep("Down-regulated", 11)),
  value = c(diff.res.num.list[["pos"]], diff.res.num.list[["neg"]]),
  res = rep(seq(0.5, 1.5, by = 0.1), 2)
)
plot.df$group <- factor(
  plot.df$group, levels = c("Up-regulated", "Down-regulated")
)

ggplot(plot.df) +
  geom_line(aes(x = res, y = value, group = group, color = group)) +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.1)) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = c(0.82, 0.88),
    legend.background = element_rect(color = "black")
  ) +
  ylab("Minimum number of DEGs across all domains") +
  xlab("Louvain resolution") +
  labs(color = "DEG type")

ggsave(
  "outputs/Visualization/fig-S7a.pdf",
  width = 5, height = 5
)

# Domain DEG (Fig. S7b) ####
## Annotate top 5 gene ####
domain.deg$label <- NA
for (domain.i in domain.list) {
  marker.idx <- which(
    domain.deg$cluster == domain.i & domain.deg$avg_log2FC > 0
  )
  anno.idx <- order(
    abs(domain.deg$avg_log2FC[marker.idx]), decreasing = TRUE
  )[1:5]
  anno.idx <- marker.idx[anno.idx]
  domain.deg[anno.idx, "label"] <- domain.deg[anno.idx, "gene"]

  marker.idx <- which(
    domain.deg$cluster == domain.i & domain.deg$avg_log2FC < 0
  )
  anno.idx <- order(
    abs(domain.deg$avg_log2FC[marker.idx]), decreasing = TRUE
  )[1:3]
  anno.idx <- marker.idx[anno.idx]
  domain.deg[anno.idx, "label"] <- domain.deg[anno.idx, "gene"]
}

ggplot(
  domain.deg,
  aes(x = cluster, y = avg_log2FC, label = label, color = cluster)
) +
  geom_jitter(position = position_jitter(seed = 1, width = 0.35)) +
  geom_text_repel(
    position = position_jitter(seed = 1, width = 0.35),
    color = "black", max.overlaps = 10
  ) +
  scale_color_manual(values = domain.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

ggsave(
  "outputs/Visualization/fig-S7b.pdf",
  width = 10, height = 5
)

# Domain DEG GO (Fig. S7c) ####
dotplot(
  deg.go,
  label_format = 100,
  font.size = 9,
  showCategory = 5
)

ggsave(
  "outputs/Visualization/fig-S7c.pdf",
  width = 12, height = 6
)

# Spatial distribution for clusters (Fig. S7d) ####
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

coo$domain <- retina.st.domain@meta.data[rownames(coo), "domain"]
pic +
  geom_point(
    data = coo,
    aes(x = coord.x, y = coord.y, fill = domain),
    size = 1, shape = 21, color = "white", stroke = 0.1
  ) +
  scale_fill_manual(values = domain.colors) +
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
  "outputs/Visualization/fig-S7d.pdf",
  width = 12, height = 4
)
