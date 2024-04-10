# Packages & settings ####
library(ggplot2)
library(ggrepel)
library(clusterProfiler)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/diff_res_num_list.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/DEG/domain_deg_go.RData")

# Domain resolution selection (Fig. S5a) ####
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
  xlab("Louvain resolution") + labs(color = "DEG type")

ggsave(
  "outputs/Visualization/fig-S5a.pdf",
  width = 5, height = 5
)

# Domain DEG (Fig. S5b) ####
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
  "outputs/Visualization/fig-S5b.pdf",
  width = 10, height = 5
)

# Domain DEG GO (Fig. S5c) ####
dotplot(
  deg.go,
  label_format = 100,
  font.size = 9,
  showCategory = 5
)

ggsave(
  "outputs/Visualization/fig-S5c.pdf",
  width = 15, height = 6
)
