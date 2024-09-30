# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)
library(viridis)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/QualityControl/retina_st_pub.RData")
load("outputs/Deconvolution/decon_results.RData")
load("outputs/Deconvolution/decon_results_pub.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")

retina.st.pub <- NormalizeData(retina.st.pub)
co.pair <- data.frame(
  gene = c("MEIS2", "PRKCA", "RRAD", "ONECUT2", "SNCG", "NRL", "HES5"),
  type = c("AC", "BC", "Cone", "HC", "RGC", "Rod", "RPC")
)[c(7, 3, 6, 5, 1, 4, 2),]

# Cell-type & marker scatter (Fig. S4a) ####
stage.mapping <- c(
  `8PCW` = 8.0, `PCW9` = 9.0, `10PCW` = 10.0, `11PCW` = 11.0,
  `PCW12_early` = 12.0, `PCW12_late` = 12.0, `13PCW` = 13.0,
  `PCW14` = 14.0, `PCW16` = 16.0, `PCW17` = 17.0
)

pic.list <- list()
for (i in seq_len(nrow(co.pair))) {
  df.1 <- data.frame(
    sample = retina.st.domain$sample.name,
    exp = retina.st.domain@assays$Spatial@data[co.pair$gene[i],],
    frac = pred.cell2loction[co.pair$type[i], Cells(retina.st.domain)],
    source = "This study"
  )
  df.2 <- data.frame(
    sample = retina.st.pub$stage,
    exp = retina.st.pub@assays$Spatial@data[co.pair$gene[i],],
    frac = pub.cell2loction[co.pair$type[i], Cells(retina.st.pub)],
    source = "Dorgau et al."
  )
  df <- rbind(df.1, df.2)
  df$sample <- stage.mapping[as.character(df$sample)]
  df <- df[df$exp > 0,]
  df <- df[df$frac > 0,]
  df <- df[sample.int(nrow(df), nrow(df)),]
  res <- summary(lm(frac ~ exp + 0, df))

  makeStars <- function(x) {
    stars <- c("****", "***", "**", "*", "ns")
    vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
    i <- findInterval(x, vec)
    stars[i]
  }


  pic.i <- ggplot(df, aes(x = exp, y = frac)) +
    geom_point(aes(color = sample, shape = source), size = 1) +
    scale_color_viridis(breaks = c(9, 11, 13, 15, 17), end = 0.9) +
    scale_shape_manual(values = c(`This study` = 15, `Dorgau et al.` = 4)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    geom_smooth(method = "lm", formula = y ~ x + 0, color = "darkblue", linewidth = 0.5) +
    ggtitle(paste0(
      "R^2=", format(res$adj.r.squared, digits = 4),
      ",", makeStars(res$coefficients[[4]])
    )) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black")
    ) +
    xlab(paste0("Expression of ", co.pair$gene[i])) +
    ylab(paste0("Fraction of ", co.pair$type[i]))

  pic.list[[i]] <- pic.i
}

(pic.list[[1]] +
  pic.list[[2]] +
  pic.list[[3]] +
  pic.list[[4]] +
  pic.list[[5]] +
  pic.list[[6]] +
  pic.list[[7]] +
  plot_spacer()) + plot_layout(ncol = 4)

ggsave(
  "outputs/Visualization/fig-S4a.pdf",
  width = 16, height = 6
)

# Cell type marker on H&E (Fig. S4b) ####
merge.coo.list <- c(
  pub.coo.list["GSM7441442_15817_8PCW_A"], coo.list["PCW9"],
  pub.coo.list["GSM7441446_15870_10PCW_A"], pub.coo.list["GSM7441450_15685_11PCW_A"],
  coo.list["PCW12_early"], coo.list["PCW12_late"],
  pub.coo.list["GSM7441454_14749_13PCW_A"], coo.list["PCW14"],
  coo.list["PCW16"], coo.list["PCW17"]
)
merge.img.list <- c(
  pub.img.list["GSM7441442_15817_8PCW_A"], img.list["PCW9"],
  pub.img.list["GSM7441446_15870_10PCW_A"], pub.img.list["GSM7441450_15685_11PCW_A"],
  img.list["PCW12_early"], img.list["PCW12_late"],
  pub.img.list["GSM7441454_14749_13PCW_A"], img.list["PCW14"],
  img.list["PCW16"], img.list["PCW17"]
)
merge.obj <- merge(retina.st.pub, retina.st.domain)

max.y <- dim(merge.img.list[[1]])[1]
start.x <- 0
for (i in 2:length(merge.img.list)) {
  start.x <- c(start.x, start.x[i - 1] + dim(merge.img.list[[i - 1]])[2] + 100)
  merge.coo.list[[i]]$coord.x <- merge.coo.list[[i]]$coord.x + start.x[i]
  if (max.y < dim(merge.img.list[[i]])[1]) {
    max.y <- dim(merge.img.list[[i]])[1]
  }
}
start.y <- NULL
for (i in seq_along(merge.coo.list)) {
  start.y <- c(start.y, ceiling((max.y - dim(merge.img.list[[i]])[1]) / 2))
  merge.coo.list[[i]]$coord.y <- merge.coo.list[[i]]$coord.y + start.y[i]
}
merge.coo <- do.call("rbind", merge.coo.list)
rownames(merge.coo) <- sub("^[^.]*\\.", "", rownames(merge.coo))

pic <- ggplot()
for (i in seq_along(merge.img.list)) {
  img.i <- merge.img.list[[i]]
  pic <- pic +
    annotation_custom(
      rasterGrob(
        img.i,
        width = unit(1, "npc"),
        height = unit(1, "npc"),
        gp = gpar(alpha = 0.8)
      ),
      xmin = start.x[i], xmax = start.x[i] + dim(img.i)[2],
      ymin = start.y[i], ymax = start.y[i] + dim(img.i)[1]
    )
}

merge.spatial.plot <- function(gene.symbol) {
  merge.coo$exp <- merge.obj@assays$Spatial@data[gene.symbol, rownames(merge.coo)]
  merge.coo$exp[merge.coo$exp > quantile(merge.coo$exp, 0.95)] <- quantile(merge.coo$exp, 0.95)
  pic.co <- pic +
    geom_point(data = merge.coo, aes(x = coord.x, y = coord.y, color = exp), size = 0.25, shape = 16) +
    scale_color_viridis(gene.symbol, option = "B") +
    coord_fixed(
      xlim = c(0, start.x[10] + dim(merge.img.list[[10]])[2]),
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

plot.list <- list()
for (gene.i in co.pair$gene) {
  plot.list[[gene.i]] <- merge.spatial.plot(gene.i)
}

plot.list[[1]] /
  plot.list[[2]] /
  plot.list[[3]] /
  plot.list[[4]] /
  plot.list[[5]] /
  plot.list[[6]] /
  plot.list[[7]]

ggsave(
  "outputs/Visualization/fig-S4b.pdf",
  width = 15, height = 10
)
