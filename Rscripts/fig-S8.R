# Packages & settings ####
library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(VennDiagram)
library(GeneOverlap)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/Domain/retina_st_domain.RData")
load("outputs/Development/reg_gene_ava.RData")
load("outputs/Development/reg_gene_df.RData")
load("outputs/DEG/domain_deg.RData")
load("outputs/DEG/reg_gene_go.RData")

# Venn of TF, RBP & DEG (Fig. S8a) ####
venn.list <- list(
  TF = tf.list,
  RBP = rbp.list,
  DEG = unique(domain.deg$gene)
)
venn.diagram(
  venn.list, disable.logging = TRUE,
  imagetype = "svg",
  fill = c("#FCB8A0", "#BDB9DA", "#7CB0D1"),
  alpha = c(0.7, 0.7, 0.7), lwd = 1,
  filename = "outputs/Visualization/fig-S8a.svg",
  width = 4, height = 4
)
print(
  intersect(
    venn.list$TF,
    intersect(venn.list$RBP, venn.list$DEG)
  )
) # HMGA1, JUN, REPIN1, ZNF207, ZNF385A, YBX1, YBX3, CXXC5, CHCHD3, SAFB 
ol.tf <- testGeneOverlap(newGeneOverlap(
  venn.list$TF, venn.list$DEG, genome.size = nrow(retina.st.domain)
))
print(paste(ol.tf@odds.ratio, ol.tf@pval))
ol.rbp <- testGeneOverlap(newGeneOverlap(
  venn.list$RBP, venn.list$DEG, genome.size = nrow(retina.st.domain)
))
print(paste(ol.rbp@odds.ratio, ol.rbp@pval))

# Number of domain associated with TF & RBP (Fig. S8b-S8c) ####
## TF ####
reg.df.tf <- reg.gene.df[reg.gene.df$type != "RBP only", ]
tf.freq <- as.data.frame(table(reg.df.tf$gene))
tf.num.freq <- as.data.frame(table(tf.freq$Freq), stringsAsFactors = FALSE)
tf.num.freq$type <- "TF"
tf.num.freq$Var1 <- factor(tf.num.freq$Var1, levels = 1:7)

pic.b <- ggplot(tf.num.freq, aes(x = Var1, y = Freq + 0.3)) +
  geom_bar(stat = "identity", fill = "#02818A") +
  geom_text(aes(label = Freq), vjust = -0.2) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("Number of associated domains") +
  ylab("Number of TFs") +
  ylim(0, 80) + coord_fixed(ratio = 0.1)

## RBP ####
reg.df.rbp <- reg.gene.df[reg.gene.df$type != "TF only", ]
rbp.freq <- as.data.frame(table(reg.df.rbp$gene))
rbp.num.freq <- as.data.frame(table(rbp.freq$Freq), stringsAsFactors = FALSE)
rbp.num.freq <- rbind(rbp.num.freq, data.frame(Var1 = "5", Freq = 0))
rbp.num.freq$type <- "RBP"
rbp.num.freq$Var1 <- factor(rbp.num.freq$Var1, levels = 1:7)

pic.c <- ggplot(rbp.num.freq, aes(x = Var1, y = Freq + 0.3)) +
  geom_bar(stat = "identity", fill ="#002672" ) +
  geom_text(aes(label = Freq), vjust = -0.2) +
  theme_classic() +
  theme(axis.text = element_text(color = "black")) +
  xlab("Number of associated domains") +
  ylab("Number of RBPs") +
  ylim(0, 80) + coord_fixed(ratio = 0.1)

pic.b | pic.c
ggsave(
  "outputs/Visualization/fig-S8bc.pdf",
  width = 8, height = 4
)

# GO network of regulatory gene (Fig. S8d) ####
cnetplot(
  reg.go, layout = "kk", showCategory = 20,
  shadowtext = "none", 
)

ggsave(
  filename = "outputs/Visualization/fig-S8d.pdf",
  width = 12, height = 12
)
