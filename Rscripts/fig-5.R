# Packages & settings ####
library(Seurat)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(GeneOverlap)
library(ggpubr)
library(ggforce)
library(ggnewscale)
library(grid)
library(viridis)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Disease/retnet.RData")
load("outputs/Disease/kegg_go_gene.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/Development/reg_gene_df.RData")
load("outputs/ImageProcessing/coordinate_HE.RData")
load("outputs/QualityControl/retina_st_qc.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Disease gene heat map (Fig. 5a) ####
exp.data <- retina.st.domain@assays$Spatial@data
domain.ds.list <- list()
for (domain.i in domain.list[1:8]) {
  domain.ds.list.i <- c()
  for (ds.i in names(retnet)) {
    ds.gene.i <- retnet[[ds.i]]
    ds.exp.i <- exp.data[ds.gene.i, retina.st.domain$domain == domain.i]
    if (length(ds.gene.i) > 1) {
      ds.exp.i <- colMeans(as.matrix(ds.exp.i))
    }
    domain.ds.list.i <- c(domain.ds.list.i, mean(ds.exp.i))
  }
  domain.ds.list[[domain.i]] <- domain.ds.list.i
}
domain.ds.df <- do.call("rbind", domain.ds.list)
colnames(domain.ds.df) <- names(retnet)

bk <- c(seq(-2, -0.01, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(
  domain.ds.df, scale = "column",
  cluster_rows = TRUE, cluster_cols = TRUE,
  treeheight_row = 30, treeheight_col = 30,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  breaks = bk, cutree_cols = 3, cutree_rows = 3,
  color = c(
    colorRampPalette(colors = c("#5D94A4", "white"))(length(bk) / 2),
    colorRampPalette(colors = c("white", "#DA3B46"))(length(bk) / 2)
  ),
  angle_col = 90, border_color = "white",
  filename = "outputs/Visualization/fig-5a.pdf",
  width = 7, height = 5
)

# Venn of Reg, DEG, & DS (Fig. 5b) ####
venn.list <- list(
  DEG = unique(domain.deg$gene),
  Reg = unique(reg.gene.df$gene),
  RetNet = unique(unlist(retnet))
)
venn.diagram(
  venn.list, disable.logging = TRUE,
  imagetype = "svg", cex = 2, cat.cex = 2,
  fill = c("#66C2A5", "#FC8D62", "#8DA0CB"),
  alpha = c(0.6, 0.6, 0.6), lwd = 3,
  filename = "outputs/Visualization/fig-5b.svg",
  width = 10, height = 10
)
print(
  intersect(
    venn.list$Reg,
    intersect(venn.list$RetNet, venn.list$DEG)
  )
)
ol.1 <- testGeneOverlap(newGeneOverlap(
  venn.list$Reg, venn.list$RetNet, genome.size = nrow(retina.st.domain)
))
print(paste(ol.1@odds.ratio, ol.1@pval))
ol.2 <- testGeneOverlap(newGeneOverlap(
  venn.list$DEG, venn.list$RetNet, genome.size = nrow(retina.st.domain)
))
print(paste(ol.2@odds.ratio, ol.2@pval))

# Retina vs non retina & expression changes (Fig. 5c-5d) ####
## Retina vs non retina ####
retina.st.qc <- NormalizeData(retina.st.qc)
exp.mat <- retina.st.qc@assays$Spatial@data

term.list <- list()
for (set.i in 1:length(key.gene.set)) {
  df.i <- data.frame(
    exp = scale(colMeans(exp.mat[key.gene.set[[set.i]],])),
    set = names(key.gene.set)[set.i],
    group = retina.st.qc$is.retina
  )
  term.list[[names(key.gene.set)[set.i]]] <- df.i

  print(wilcox.test(
    x = df.i[df.i$group == "Retina", "exp"],
    y = df.i[df.i$group == "Non retina", "exp"],
    alternative = "greater"
  )$p.value)
}
key.gene.df <- do.call("rbind", term.list)

pic.c <- ggplot(
  key.gene.df,
  aes(x = set, y = exp, fill = group)
) +
  geom_boxplot(
    outlier.alpha = 0, size = 0.1
  ) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", method = "wilcox.test",
    method.args = list(alternative = "greater")
  ) +
  scale_y_continuous(limits = c(-2.5, 3.5)) +
  scale_fill_manual(values = c(
    `Retina` = "#A6D854", `Non retina` = "#FFD92F"
  )) +
  theme_bw() +
  labs(fill = "Tissue") +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Normalized avg. expr.")

## Gene set expression changes ####
exp.mat <- retina.st.domain@assays$Spatial@data

term.list <- list()
for (set.i in 1:length(key.gene.set)) {
  df.i <- data.frame(
    exp = scale(colMeans(exp.mat[key.gene.set[[set.i]],])),
    set = names(key.gene.set)[set.i],
    sample = retina.st.domain$sample.name
  )
  term.list[[names(key.gene.set)[set.i]]] <- df.i
}
key.gene.df <- do.call("rbind", term.list)

pic.d <- ggplot(
  key.gene.df,
  aes(x = sample, y = exp, group = set, color = set)
) +
  geom_smooth(
    method = "loess", formula = y ~ x,
    span = 1, se = FALSE, size = 0.5
  ) +
  scale_color_brewer(palette = "Set3") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  ylab("Normalized avg. expr.") +
  labs(color = "Gene set")

pic.c | pic.d
ggsave(
  "outputs/Visualization/fig-5cd.pdf",
  width = 7, height = 4
)

# VAM spatial expression (Fig. 5e) ####
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
coo$exp <- colMeans(retina.st.domain@assays$Spatial@data[key.gene.set$VAM, rownames(coo)])
coo$exp[coo$exp > quantile(coo$exp, 0.95)] <- quantile(coo$exp, 0.95)

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

pic +
  geom_point(
    data = coo,
    aes(x = coord.x, y = coord.y, fill = exp),
    size = 0.9, shape = 21, stroke = 0.1, color = "white"
  ) +
  scale_fill_viridis("VAM", option = "B") +
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
    panel.spacing = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "outputs/Visualization/fig-5e.pdf",
  width = 10, height = 5
)
