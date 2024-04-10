# Packages & settings ####
library(Seurat)
library(ggplot2)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Development/reg_gene_ava.RData")
load("outputs/Development/reg_gene_df.RData")

retina.st.domain$domain <- factor(
  domain.list[as.character(retina.st.domain$domain)], levels = domain.list
)
reg.gene.df$domain <- factor(reg.gene.df$domain, levels = domain.list)
reg.gene.df$source <- factor(reg.gene.df$source, levels = sample.list)
reg.gene.df$target <- factor(reg.gene.df$target, levels = sample.list)

# Top key TF/RBP visualization (Fig. 4a) ####
## Select top TF/RBP ####
top.k <- 5
reg.unique <- unique(reg.gene.df[, c("domain", "source")])
reg.sub.list <- list()
for (reg.i in 1:nrow(reg.unique)) {
  reg.sub.idx.i <- (
    (reg.gene.df$domain == reg.unique$domain[reg.i]) &
      (reg.gene.df$source == reg.unique$source[reg.i])
  )
  reg.sub.i <- reg.gene.df[reg.sub.idx.i, ]
  
  reg.sub.top.i <- order(
    reg.sub.i$avg_log2FC, decreasing = TRUE
  )[1:min(top.k, nrow(reg.sub.i))]
  reg.sub.i <- reg.sub.i[reg.sub.top.i, ]
  reg.sub.list <- c(reg.sub.list, list(reg.sub.i))
}
reg.sub.df <- do.call("rbind", reg.sub.list)

## Re-order ####
reg.sub.df.order <- order(
  reg.sub.df$source, reg.sub.df$target, reg.sub.df$domain
)
reg.sub.df <- reg.sub.df[reg.sub.df.order, ]

reg.sub.df$name <- paste(
  reg.sub.df$source, reg.sub.df$target, reg.sub.df$domain
)
reg.sub.df$name <- factor(
  reg.sub.df$name,
  levels = unique(reg.sub.df$name)[length(unique(reg.sub.df$name)):1]
)

reg.sub.gene.unique <- unique(reg.sub.df$gene)
reg.sub.df$gene <- factor(
  reg.sub.df$gene,
  levels = reg.sub.gene.unique
)

## Refined dot plot ####
gene.color <- unique(data.frame(
  gene = reg.sub.df$gene, type = reg.sub.df$type, color = "NA"
))
rownames(gene.color) <- gene.color$gene
gene.color <- gene.color[reg.sub.gene.unique, ]

gene.color$color[gene.color$type == "TF only"] <- "#02818A"
gene.color$color[gene.color$type == "RBP only"] <- "#002672"
gene.color$color[gene.color$type == "Both"] <- "#6DEEEE"

ggplot(
  reg.sub.df,
  aes(x = gene, y = name, color = avg_log2FC, size = -log10(p_val_adj))
) +
  geom_point() +
  scale_color_gradient(
    low = "#BFD3E6", high = "#88419D"
  ) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.border = element_rect(
      fill = NA, color = "black", size = 1, linetype = "solid"
    ),
    axis.text.x = element_text(
      color = gene.color$color, angle = 90, hjust = 0,
      face = "bold"
    ),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) + coord_fixed() +
  labs(size = "-log10(p.adj)", color = "avg_log2FC")

ggsave(
  "outputs/Visualization/fig-4a.pdf",
  width = 12, height = 8
)

# Key TF/RBP example (Fig. 4b) ####
## Group gene expression ####
key.gene <- c("LGALS1", "DCN", "NRL", "NEUROD1", "LDHA", "HES5")
term.list <- list()
exp.mat <- retina.st.domain@assays$Spatial@data[key.gene, ]
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

key.df$label <- as.character(key.df$domain)
key.df$label[key.df$sample != "PCW17"] <- ""
key.df$label.loc <- as.integer(key.df$domain) - 0.5

## Create donut plot ####
max.num <- 2
min.num <- -2
pic.list <- list()
for (gene.i in key.gene) {
  key.df.i <- key.df[key.df$gene == gene.i, ]
  key.df.i$values <- scale(key.df.i$values)
  key.df.i$values[key.df.i$values > max.num] <- max.num
  key.df.i$values[key.df.i$values < min.num] <- min.num
  
  key.df.blank.i <- key.df.i[key.df.i$sample == "PCW9", ]
  key.df.blank.i$sample <- "blank"
  key.df.blank.i$values <- 0
  key.df.i <- rbind(key.df.i, key.df.blank.i)
  key.df.i$sample <- factor(
    key.df.i$sample, levels = c("blank", sample.list)
  )
  
  pic.list[[gene.i]] <- ggplot(key.df.i) +
    geom_bar(
      aes(x = sample, y = gene, fill = values),
      stat = "identity", size = 1, color = "white"
    ) +
    scale_fill_gradient2(
      low = "#5D94A4", mid = "white", high = "#DA3B46",
      breaks = c(-2, -1, 0, 1, 2), limits = c(min.num, max.num)
    ) +
    geom_text(
      aes(x = sample, y = label.loc, label = label),
      color = "black"
    ) +
    coord_polar(theta = "y", direction = -1) +
    theme_void() + ggtitle(gene.i)
}

(pic.list[["LGALS1"]] | pic.list[["DCN"]] | pic.list[["NRL"]]) /
  (pic.list[["NEUROD1"]] | pic.list[["LDHA"]] | pic.list[["HES5"]])
ggsave(
  "outputs/Visualization/fig-4b.pdf",
  width = 12, height = 8
)

# VEGFA/PGF visualization (Fig. 4c) ####
key.gene <- c("VEGFA", "PGF")
term.list <- list()
exp.mat <- retina.st.domain@assays$Spatial@data[key.gene, ]
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
    geom_smooth(method = "loess", formula = y ~ x, span = 2, se = FALSE) +
    scale_color_manual(values = domain.colors) +
    scale_y_continuous(
      na.value = 0,
      limits = c(max(0, min(df$values - 1e-2)), max(df$values + 1e-2))
    ) +
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

pic.c1 <- key.gene.plot(key.df[key.df$gene == "VEGFA", ]) +
  theme(legend.position = "none") +
  ylab("Average expression")
pic.c2 <- key.gene.plot(key.df[key.df$gene == "PGF", ]) +
  theme(axis.title.y = element_blank())

pic.c1 | pic.c2
ggsave(
  "outputs/Visualization/fig-4c.pdf",
  width = 6, height = 5
)
