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
load("outputs/Disease/kegg_go_gene.RData")

# Spatial plot function ####
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

spatial.gene.set.plot <- function(gene.set, name) {
  coo$exp <- colMeans(retina.st.domain@assays$Spatial@data[gene.set, rownames(coo)])
  coo$exp[coo$exp > quantile(coo$exp, 0.95)] <- quantile(coo$exp, 0.95)
  pic.i <- pic +
    geom_point(
      data = coo,
      aes(x = coord.x, y = coord.y, fill = exp),
      size = 0.7, shape = 21, stroke = 0.05, color = "white"
    ) +
    scale_fill_viridis(name, option = "B") +
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
  pic.i
}

# Gene sets spatial expression (Fig. 14) ####
pic.a <- spatial.gene.set.plot(key.gene.set$FAM, "FAM")
pic.b <- spatial.gene.set.plot(key.gene.set$FSVB, "FSVB")
pic.c <- spatial.gene.set.plot(key.gene.set$FSVC, "FSVC")
pic.d <- spatial.gene.set.plot(key.gene.set$VDM, "VDM")
pic.e <- spatial.gene.set.plot(key.gene.set$VEM, "VEM")
pic.f <- spatial.gene.set.plot(key.gene.set$VKM, "VKM")

pic.a / pic.b / pic.c / pic.d / pic.e / pic.f

ggsave(
  "outputs/Visualization/fig-S14.pdf",
  width = 10, height = 12
)
