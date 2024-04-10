# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggnewscale)
library(patchwork)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Deconvolution/decon_results.RData")

# Cell-type & marker co-localization (Fig. S3) ####
umap.co.plot <- function(data) {
  plt <- ggplot(data) +
    geom_point(
      aes(x = UMAP_1, y = UMAP_2, color = `Exp.`),
      size = 2
    ) +
    scale_color_gradient(
      low = "#D6D6D6", high = "#EE692E",
      guide = guide_colorbar(order = 2)
    ) +
    new_scale_color() +
    geom_point(
      aes(x = UMAP_1, y = UMAP_2, color = `Frac.`),
      size = 0.2
    ) +
    scale_color_gradient(
      low = "#D6D6D6", high = "#9A1DEE",
      guide = guide_colorbar(order = 1)
    ) +
    coord_fixed(ratio = 2) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank()
    ) +
    theme(
      legend.position = "right",
      legend.key.width = unit(0.4, "cm"),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    )
  plt
}

co.pair <- data.frame(
  gene = c("MEIS2", "PRKCA", "RRAD", "ONECUT2", "SNCG", "NRL", "HES5"),
  type = c("AC", "BC", "Cone", "HC", "RGC", "Rod", "RPC")
)[c(7, 3, 6, 5, 1, 4, 2), ]

plot.list <- list()
df <- data.frame(
  UMAP_1 = retina.st.domain@reductions$umap@cell.embeddings[, 1],
  UMAP_2 = retina.st.domain@reductions$umap@cell.embeddings[, 2]
)
for (i in 1:nrow(co.pair)) {
  df$`Exp.` <- retina.st.domain@assays$Spatial@data[co.pair$gene[i], ]
  df$`Frac.` <- pred.cell2loction[co.pair$type[i], Cells(retina.st.domain)]
  
  plot.list[[i]] <- umap.co.plot(df)
}

(plot.list[[1]] | plot.list[[2]] | plot.list[[3]] | plot.list[[4]]) /
  (plot.list[[5]] | plot.list[[6]] | plot.list[[7]] | plot_spacer())

ggsave(
  filename = "outputs/Visualization/fig-S3.pdf",
  width = 20, height = 10
)
