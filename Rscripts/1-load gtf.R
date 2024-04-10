# Packages & settings ####
rm(list = ls())
source("scripts/utils.R")

# Split GTF file ####
ori.ref.gtf.file <- "data/Reference/GRCh38-GENCODE_v43/genes.gtf.gz"
ref.gtf <- read.csv(
  file = ori.ref.gtf.file,
  header = FALSE,
  sep = "\t",
  skip = 5
)

ref.gtf <- ref.gtf[ref.gtf$V3 == "gene", ]
ref.gtf.info <- strsplit(ref.gtf$V9, split = "[ ;]")
ref.gtf[["gene.id"]] <- sapply(X = ref.gtf.info, FUN = "[", 2)
ref.gtf[["gene.type"]] <- sapply(X = ref.gtf.info, FUN = "[", 8)
ref.gtf[["gene.symbol"]] <- sapply(X = ref.gtf.info, FUN = "[", 11)
rownames(ref.gtf) <- ref.gtf$gene.id

# Convert duplicated gene symbols ####
ref.gtf <- ref.gtf[, c(1, 4, 5, 7, 10, 12, 11)]
colnames(ref.gtf) <- c(
  "chromosome", "start", "end", "strand",
  "gene.id", "gene.symbol", "gene.type"
)

duplicated.idx <- which(duplicated(ref.gtf$gene.symbol))
num.loop <- 1
while (length(duplicated.idx) > 0) {
  old.symbol <- ref.gtf$gene.symbol[duplicated.idx]
  if (num.loop > 1) {
    old.symbol <- sapply(strsplit(old.symbol, "[.]"), "[", 1)
  }
  new.symbol <- paste0(old.symbol, ".", num.loop)
  
  ref.gtf$gene.symbol[duplicated.idx] <- new.symbol
  duplicated.idx <- which(duplicated(ref.gtf$gene.symbol))
  num.loop <- num.loop + 1
}

# Screening out genes on chromosome Y ####
chromosome.ava <- paste0("chr", c(1:22, "M", "X"))
gene.idx <- ref.gtf$chromosome %in% chromosome.ava
ref.gtf <- ref.gtf[gene.idx, ]

# Save ####
write.table(
  ref.gtf,
  file = ref.gtf.file,
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
