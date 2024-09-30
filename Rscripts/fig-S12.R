# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)

# Growth factor visualization (Fig. S12) ####
key.gene <- c(
  "IGF1", "IGF1R", "IGF2R",
  "FGF1", "FGFR1", "FGFR2",
  "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
  "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2",
  "NGF", "BDNF", "NTF3",
  "EGF", "EGFR"
)
term.list <- list()
exp.mat <- retina.st.domain@assays$Spatial@data[key.gene,]
for (gene.i in key.gene) {
  for (domain.i in domain.list[1:8]) {
    for (sample.i in sample.list) {
      idx.i <- (
        (retina.st.domain$sample.name == sample.i) &
          (retina.st.domain$domain == domain.i)
      )
      term.list <- c(
        term.list,
        list(data.frame(
          values = mean(exp.mat[gene.i, idx.i]),
          domain = domain.i,
          sample = sample.i,
          gene = gene.i
        ))
      )
    }
  }
}
key.df <- do.call("rbind", term.list)
key.df$values[is.na(key.df$values)] <- 0

key.df$domain <- factor(key.df$domain, levels = domain.list[1:8])
key.df$sample <- factor(key.df$sample, levels = sample.list)

key.gene.plot <- function(df) {
  pic <- ggplot(
    df,
    aes(x = sample, y = values, color = domain, group = domain)
  ) +
    geom_smooth(
      method = "loess", formula = y ~ x, span = 2,
      se = FALSE, linewidth = 0.5
    ) +
    scale_color_manual(values = domain.colors) +
    scale_y_continuous(na.value = 0) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    facet_grid(. ~ gene)
  pic
}

pic.list <- list()
for (gene.i in unique(key.df$gene)) {
  pic.list <- c(pic.list, list(
    key.gene.plot(key.df[key.df$gene == gene.i,]) +
      theme(legend.position = "none") +
      ylab("Average expression")
  ))
}

(pic.list[[1]] |
  pic.list[[2]] |
  pic.list[[3]] |
  pic.list[[4]]) /
  (pic.list[[5]] |
    pic.list[[6]] |
    pic.list[[7]] |
    pic.list[[8]]) /
  (pic.list[[9]] |
    pic.list[[10]] |
    pic.list[[11]] |
    pic.list[[12]]) /
  (pic.list[[13]] |
    pic.list[[14]] |
    pic.list[[15]] |
    pic.list[[16]]) /
  (pic.list[[17]] |
    pic.list[[18]] |
    pic.list[[19]] |
    pic.list[[20]])

ggsave(
  "outputs/Visualization/fig-S12.pdf",
  width = 12, height = 12
)
