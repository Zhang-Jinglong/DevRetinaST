# Packages & settings ####
library(ggplot2)
library(ggpubr)
library(igraph)
library(ggraph)

rm(list = ls())
source("scripts/utils.R")

# SCRIPro strength distribution (Fig. S11a) ####
tf.target.net <- read.csv(
  "outputs/Development/SCRIPro/scripro_tf_target_strength.csv"
)[, 2:4]

strength.idx.1 <- tf.target.net$Strength > 0.3
tf.target.net.sub.1 <- tf.target.net[strength.idx.1,]
tf.freq.1 <- as.data.frame(table(tf.target.net.sub.1$TF))
tf.freq.1$range <- "x > 0.3"
strength.idx.2 <- (tf.target.net$Strength > 0.25) & (tf.target.net$Strength <= 0.3)
tf.target.net.sub.2 <- tf.target.net[strength.idx.2,]
tf.freq.2 <- as.data.frame(table(tf.target.net.sub.2$TF))
tf.freq.2$range <- "0.25 < x <= 0.3"
strength.idx.3 <- (tf.target.net$Strength > 0.2) & (tf.target.net$Strength <= 0.25)
tf.target.net.sub.3 <- tf.target.net[strength.idx.3,]
tf.freq.3 <- as.data.frame(table(tf.target.net.sub.3$TF))
tf.freq.3$range <- "0.2 < x <= 0.25"
tf.freq <- do.call("rbind", list(tf.freq.1, tf.freq.2, tf.freq.3))
tf.freq$Var1 <- as.character(tf.freq$Var1)
tf.order <- tf.freq$Var1[order(tf.freq$range, tf.freq$Freq, decreasing = TRUE)]
tf.order <- tf.order[!duplicated(tf.order)]
tf.sum <- aggregate(tf.freq$Freq, by = list(tf.freq$Var1), sum)
tf.order <- tf.sum$Group.1[order(tf.sum$x, decreasing = TRUE)]
tf.freq$Gene <- factor(tf.freq$Var1, levels = tf.order)

ggplot(tf.freq) +
  geom_bar(
    aes(x = Gene, y = Freq, fill = range),
    position = "stack", stat = "identity"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "outputs/Visualization/fig-S11a.pdf",
  width = 12, height = 3
)

# GRN (Fig. S11b) ####
tf.list <- c("NR1D1", "EPAS1", "HES4", "POU4F1", "CRX")
tf.target.net.sub.1 <- tf.target.net.sub.1[tf.target.net.sub.1$TF %in% tf.list,]
nodes <- data.frame(
  Gene = c(tf.list, unique(tf.target.net.sub.1$Target)),
  Type = c(rep("TF", length(tf.list)), rep("Regulon", length(unique(tf.target.net.sub.1$Target))))
)
nodes <- nodes[!duplicated(nodes$Gene),]
rownames(nodes) <- nodes$Gene
graph <- graph_from_data_frame(tf.target.net.sub.1, directed = TRUE, vertices = nodes)

lay <- create_layout(graph, layout = "fr")
ggraph(lay) +
  geom_edge_link(
    aes(edge_width = Strength),
    color = "grey",
    arrow = arrow(type = "closed", length = unit(2, "mm")),
    start_cap = circle(3, "mm"),
    end_cap = circle(3, "mm")
  ) +
  geom_edge_loop(
    aes(edge_width = Strength),
    color = "grey",
    arrow = arrow(type = "closed", length = unit(2, "mm")),
    start_cap = circle(3, "mm"),
    end_cap = circle(3, "mm")
  ) +
  geom_node_point(
    aes(color = Type),
    size = 5
  ) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_brewer(palette = "Set1") +
  scale_edge_width(range = c(0.5, 1.5), breaks = c(0.32, 0.36, 0.40, 0.42, 0.44)) +
  theme_void()

ggsave(
  "outputs/Visualization/fig-S11b.pdf",
  width = 8, height = 8
)
