# Packages & settings ####
library(Seurat)
library(ggplot2)
library(grid)
library(viridis)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")
load("outputs/MergeData/retina_bulk_pub.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Growth factor spatial expression (Fig. S13a) ####
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

spatial.plot <- function(gene.symbol) {
  coo$exp <- retina.st.domain@assays$Spatial@data[gene.symbol, rownames(coo)]
  coo$exp[coo$exp > quantile(coo$exp, 0.95)] <- quantile(coo$exp, 0.95)
  pic.co <- pic +
    geom_point(data = coo, aes(x = coord.x, y = coord.y, color = exp), size = 0.35, shape = 16) +
    scale_color_viridis(option = "B") +
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
  pic.co
}

key.gene <- c("VEGFA", "PGF", "CCN2", "IGF2", "FGF2", "VGF")
plot.list <- list()
for (gene.i in key.gene) {
  plot.list[[gene.i]] <- spatial.plot(gene.i)
}

plot.list[[1]] /
  plot.list[[2]] /
  plot.list[[3]] /
  plot.list[[4]] /
  plot.list[[5]] /
  plot.list[[6]]

ggsave(
  "outputs/Visualization/fig-S13a.pdf",
  width = 10, height = 10
)

# Growth factor bulk organoid expression (Fig. S13b) ####
key.gene <- c("VEGFA", "PGF", "CCN2", "IGF2", "FGF2", "VGF")
retina.bulk <- retina.bulk[key.gene,]
retina.bulk$gene <- rownames(retina.bulk)
retina.bulk <- reshape2::melt(retina.bulk, id.vars = "gene")
retina.bulk$stage <- as.integer(substr(
  sapply(strsplit(as.character(retina.bulk$variable), "_"), "[", 1),
  start = 2, stop = 4
))
retina.bulk <- retina.bulk[retina.bulk$stage != 200,]
retina.bulk$gene <- factor(retina.bulk$gene, levels = key.gene)

ggplot(
  retina.bulk,
  aes(x = stage, y = value)
) +
  geom_point(size = 1) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 2,
    se = TRUE, linewidth = 0.5
  ) +
  facet_wrap(. ~ gene, nrow = 1) +
  scale_x_continuous(limits = c(40, 130), breaks = c(40, 60, 80, 100, 120)) +
  annotate(
    "rect",
    xmin = 55, xmax = 105, ymin = 1, ymax = max(retina.bulk$value) + 1000,
    fill = "darkred", alpha = 0.1
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black")
  )

ggsave(
  "outputs/Visualization/fig-S13b.pdf",
  width = 12, height = 3
)
