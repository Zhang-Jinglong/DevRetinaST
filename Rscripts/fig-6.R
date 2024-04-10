# Packages & settings ####
library(Seurat)
library(ggplot2)
library(ggrepel)
library(ggplot2)
library(ggnewscale)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/ImageProcessing/coordinate.RData")

spot.list <- Cells(retina.st.domain)
retina.st.domain$coord.x <- spatial.coord.new[spot.list, "coord.x"]
retina.st.domain$coord.y <- spatial.coord.new[spot.list, "coord.y"]

lr.table <- read.csv(
  "outputs/CCC/COMMOT/retina_db_filter.csv", row.names = 1
)
rownames(lr.table) <- paste(lr.table$X0, lr.table$X1, sep = "-")

lr.summary.list <- list()
for (sample.i in sample.list) {
  lr.summary.i <- read.csv(
    paste0("outputs/CCC/COMMOT/commot_", sample.i, "_lr.csv"),
    header = FALSE
  )
  colnames(lr.summary.i) <- c("LR", "value")
  lr.summary.i$sample <- sample.i
  lr.summary.i$value <- as.numeric(scale(lr.summary.i$value))
  
  lr.summary.list[[sample.i]] <- lr.summary.i
}
lr.summary <- do.call("rbind", lr.summary.list)
lr.summary.df <- aggregate(
  lr.summary$value, by = list(LR = lr.summary$LR), mean
)
lr.summary.df$type <- lr.table[lr.summary.df$LR, "X3"]

# Spatial communication score (Fig. 6b) ####
lr.summary.df <- lr.summary.df[order(lr.summary.df$type, lr.summary.df$LR), ]
lr.summary.df$LR <- factor(lr.summary.df$LR, levels = unique(lr.summary.df$LR))
lr.summary.df$label <- as.character(lr.summary.df$LR)
lr.summary.df$label[lr.summary.df$x < 2.0] <- ""

ggplot(lr.summary.df) +
  geom_point(aes(x = LR, y = x, color = type)) +
  scale_color_manual(values = c(
    `Cell-Cell Contact` = "#658A4B", `Secreted Signaling` = "#53379C"
  )) +
  geom_text_repel(
    aes(x = LR, y = x, label = label), color = "black"
  ) +
  geom_hline(yintercept = 2.0, linetype = 3, color = "red") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  xlab("Ligand-receptor pair") + labs(color = "LR type") +
  ylab("Average scaled SCS")

ggsave(
  "outputs/Visualization/fig-6b.pdf",
  width = 7.5, height = 5
)

# MDK & LRP1 co-localization (Fig. 6c) ####
plot.data <- data.frame(
  l.value = retina.st.domain@assays$Spatial@data["MDK", ],
  r.value = retina.st.domain@assays$Spatial@data["LRP1", ],
  x = (
    retina.st.domain$coord.x +
      coord.add.2[retina.st.domain$sample.name, "x"]
  ),
  y = (
    retina.st.domain$coord.y +
      coord.add.2[retina.st.domain$sample.name, "y"]
  )
)
plot.data$l.value <- scale(plot.data$l.value)
plot.data$l.value[plot.data$l.value > 2.5] <- 2.5
plot.data$l.value[plot.data$l.value < -2.5] <- -2.5
plot.data$r.value <- scale(plot.data$r.value)
plot.data$r.value[plot.data$r.value > 2.5] <- 2.5
plot.data$r.value[plot.data$r.value < -2.5] <- -2.5

pic.c <- ggplot()
for (sample.i in sample.list) {
  pic.c <- create.empty.bk(
    pic.c, coord.add.2[sample.i, "x"], coord.add.2[sample.i, "y"]
  )
}

pic.c +
  geom_point(
    data = plot.data, size = 0.5,
    mapping = aes(x = x, y = y, color = l.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#EE692E",
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_color() +
  geom_point(
    data = plot.data, size = 0.1,
    mapping = aes(x = x, y = y, color = r.value)
  ) +
  scale_color_gradient2(
    low = "#8C9AA4", mid = "#CCCCCC", high = "#9A1DEE",
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_color() +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  coord_fixed()

ggsave(
  "outputs/Visualization/fig-6c.pdf",
  width = 6, height = 5
)
