# Packages & settings ####
library(Seurat)
library(WGCNA)
library(preprocessCore)
options(stringsAsFactors = FALSE)

rm(list = ls())
source("scripts/utils.R")

# Load data ####
load("outputs/QualityControl/retina_st_filter.RData")

# Function to construct the network and identify modules
construct.network <- function(
    datExprdataOne, power, networkType,
    corType, output.dir
) {
  blockwiseModules(
    datExprdataOne, power = power,
    TOMType = networkType, minModuleSize = 20, deepSplit = 3,
    trapErrors = TRUE,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE, corType = corType,
    maxPOutliers = 0.1, loadTOMs = TRUE,
    saveTOMFileBase = paste0(output.dir, "tom"),
    verbose = 5
  )
}

# Define a function to perform WGCNA normalization and analysis
wgcna.analysis <- function(wgcna.data, threshold.p, output.dir) {
  # Enable multiple-threading in WGCNA
  allowWGCNAThreads()
  
  # Normalize the data using quantile normalization
  wgcna.normalize <- normalize.quantiles(wgcna.data)
  rownames(wgcna.normalize) <- rownames(wgcna.data)
  colnames(wgcna.normalize) <- colnames(wgcna.data)
  
  # Transpose the normalized data for WGCNA analysis
  datExprdataOne <- t(wgcna.normalize)
  
  # Filter out the genes with zero expression in more than threshold.p
  threshold.s <- threshold.p * nrow(datExprdataOne)
  gene.filter <- colSums(datExprdataOne < 1e-7) > threshold.s
  datExprdataOne <- datExprdataOne[, !gene.filter]
  
  # Check for the presence of good samples and genes
  gsg <- goodSamplesGenes(datExprdataOne, verbose = 3)
  if (!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed
    if (sum(!gsg$goodGenes) > 0) {
      printFlush(paste(
        "Removing genes:",
        paste(colnames(datExprdataOne)[!gsg$goodGenes], collapse = ",")
      ))
    }
    if (sum(!gsg$goodSamples) > 0) {
      printFlush(paste(
        "Removing samples:",
        paste(rownames(datExprdataOne)[!gsg$goodSamples], collapse = ",")
      ))
    }
    # Remove the offending genes and samples from the data
    datExprdataOne <- datExprdataOne[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Define network parameters
  networkType <- "signed"
  corType <- "bicor"
  
  # Choose the power for the soft-threshold
  powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
  sft <- pickSoftThreshold(
    datExprdataOne, powerVector = powers,
    networkType = networkType, verbose = 2
  )
  save(sft, file = paste0(output.dir, "sft.RData"))
  
  # Select the power value
  power <- sft$powerEstimate
  
  # Construct the network and identify modules  
  net <- construct.network(
    datExprdataOne, power, networkType, corType, output.dir
  )
  
  # Save eigenvalue for each modules #
  moduleLabels <- as.data.frame(net$colors)
  moduleLabels <- cbind(moduleLabels, rownames(moduleLabels))
  module_number <- length(unique(net$colors))
  for (i in c(1:(module_number - 1))) {
    gene <- moduleLabels[which(moduleLabels == i), 2]
    write.table(
      gene,
      paste0(output.dir, "ME", as.character(i), ".txt"),
      row.names = FALSE, quote = FALSE, col.names = FALSE
    )
  }
  
  # Save the results
  save(net, file = paste0(output.dir, "net.RData"))
}

retina.st.filter <- NormalizeData(retina.st.filter)
wgcna.data <- as.matrix(retina.st.filter@assays$Spatial@data)

set.seed(3407)
wgcna.analysis(wgcna.data, 0.8, "outputs/WGCNA/")
