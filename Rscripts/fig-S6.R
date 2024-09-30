# Packages & settings ####
library(Seurat)
library(ggplot2)
library(patchwork)

rm(list = ls())
source("scripts/utils.R")

# UMAP visualization of embeddings from 5 methods (Fig. S6a) ####
load("outputs/Domain/retina_st_domain.RData")
retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2],
  Harmony = retina.st.domain$domain
)
pic.harmony <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = Harmony
  ), size = 0.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_manual(values = domain.colors) +
  coord_fixed(ratio = 2)

load("outputs/Domain/retina_st_domain_bass.RData")
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap.bass@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap.bass@cell.embeddings[, 2],
  BASS = paste0("D", as.character(retina.st.domain$domain))
)
pic.bass <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = BASS
  ), size = 0.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_brewer(palette = "Set3") +
  coord_fixed(ratio = 1.5)

load("outputs/Domain/retina_st_domain_graphst.RData")
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap.graphst@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap.graphst@cell.embeddings[, 2],
  GraphST = paste0("D", as.character(retina.st.domain$domain))
)
pic.graphst <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = GraphST
  ), size = 0.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_brewer(palette = "Set3") +
  coord_fixed(ratio = 1)

load("outputs/Domain/retina_st_domain_cca.RData")
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap.cca@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap.cca@cell.embeddings[, 2],
  Seurat_CCA = paste0("D", as.integer(retina.st.domain$domain))
)
pic.cca <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = Seurat_CCA
  ), size = 0.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_brewer("Seurat-CCA", palette = "Set3") +
  coord_fixed(ratio = 1)

load("outputs/Domain/retina_st_domain_scvi.RData")
umap.data <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap.scvi@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap.scvi@cell.embeddings[, 2],
  scVI = paste0("D", as.integer(retina.st.domain$domain))
)
pic.scvi <- ggplot(umap.data) +
  geom_point(aes(
    x = UMAP_1, y = UMAP_2,
    color = scVI
  ), size = 0.5) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  scale_color_brewer(palette = "Set3") +
  coord_fixed(ratio = 1)

(pic.harmony |
  pic.bass |
  pic.graphst |
  pic.cca |
  pic.scvi) +
  plot_layout(widths = c(2, 2, 2, 2, 2))

ggsave(
  filename = "outputs/Visualization/fig-S6a.pdf",
  width = 15, height = 4
)

# iLISI & accuracy (Fig. S6b-S6c) ####
load("outputs/Domain/lisi_6_methods.RData")
lisi$method <- factor(lisi$method, levels = c("Raw", "Harmony", "BASS", "GraphST", "Seurat-CCA", "scVI"))
pic.lisi <- ggplot(lisi) +
  geom_boxplot(
    aes(x = method, y = iLISI, fill = method),
    outlier.size = 1
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  scale_fill_brewer(palette = "Set2")

acc <- read.csv("outputs/Domain/CrossVal/acc_score.csv")
acc$method <- factor(acc$method, levels = c("Harmony", "BASS", "GraphST", "Seurat-CCA", "scVI"))
pic.acc <- ggplot(acc) +
  geom_boxplot(
    aes(y = Accuracy, fill = Classifier),
    outlier.size = 1
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  scale_fill_brewer(palette = "Accent") +
  facet_grid(~method)

pic.lisi | pic.acc

ggsave(
  "outputs/Visualization/fig-S6bc.pdf",
  width = 13.5, height = 5
)
