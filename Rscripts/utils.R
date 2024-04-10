# Reference ####
ref.gtf.file <- "outputs/Reference/genes.gtf.csv"

# Spot category list ####
slice.list <- c("PCW9_12", "PCW12", "PCW14", "PCW16", "PCW17")
sample.list <- c("PCW9", "PCW12_early", "PCW12_late", "PCW14", "PCW16", "PCW17")
radial.list <- c(
  paste0("NBL ", 1:4, "/4"), paste0("GCL ", 1:2, "/2")
)
cell.type.list <- c(
  "Cone", "Rod", "AC", "BC", "HC", "RGC",
  "RPC", "Neurogenic", "BC-Photo Precur", "AC-HC Precur"
)
domain.list <- c(
  `0` = "D1", `2` = "D2", `6` = "D3", `5` = "D4", `1` = "D5",
  `3` = "D6", `4` = "D7", `7` = "D8", `8` = "D9"
)

# Color palette ####
sample.colors <- c(
  `PCW9` = "#E2DE2C", `PCW12_early` = "#EEBAA2", `PCW12_late` = "#FF7F00",
  `PCW14` = "#A65629", `PCW16` = "#F781BF", `PCW17` = "#984EA3"
)
radial.colors <- c(
  `NBL 1/4` = "#E31A1C", `NBL 2/4` = "#FF5560",
  `NBL 3/4` = "#FB91D6", `NBL 4/4` = "#CAB2D6",
  `GCL 1/2` = "#A6CEE3", `GCL 2/2` = "#1F78B4"
)
cell.type.colors <- c(
  `Cone` = "#3066B8", `Rod` = "#4DAF4A",
  `AC` = "#FF7F00", `BC` = "#6A3D9A", `HC` = "#9A2202", `RGC` = "#F0A4AE",
  `RPC` = "#A1004E", `Neurogenic` = "#D765AC",
  `BC-Photo Precur` = "#22818A", `AC-HC Precur` = "#BF5B16"
)
domain.colors <- c(
  `D1` = "#A9170E", `D2` = "#E25B5F", `D3` = "#FEB7BB",
  `D4` = "#DBCAE4", `D5` = "#906C9C", `D6` = "#74B5D6",
  `D7` = "#718CB2", `D8` = "#6862A7", `D9` = "#3F3C64"
)
module.colors <- c(
  `ME1` = "#3AA887", `ME2` = "#BF5B17", `ME3` = "#9490C0", `ME4` = "#E6AC02"
)

# Other constants ####
coord.add.1 <- data.frame(
  x = seq(0, 5, 1) * 17,
  y = rep(0, 6),
  row.names = sample.list
)
coord.add.2 <- data.frame(
  x = rep(seq(0, 2, 1) * 17, 2),
  y = c(rep(0, 3), rep(-16, 3)),
  row.names = sample.list
)

# Create spatial graph ####
create.spatial.graph <- function(coordinate) {
  distance <- as.matrix(
    dist(coordinate[, 4:5], method = "euclidean")
  ) + (diag(nrow(coordinate)) * 1e7)
  
  distance[distance > min(distance) * 1.5] <- 0
  distance[distance > 0] <- 1
  
  distance
}

# Create label mesh ####
create.label.mesh <- function(labels, start) {
  point.old <- start
  labels$flag <- FALSE
  count <- 0
  
  labels.mesh <- list()
  while (count < (nrow(labels))) {
    dist.i <- sqrt(
      ((labels$x - point.old$x) ^ 2) + ((labels$y - point.old$y) ^ 2)
    )
    dist.i[labels$flag] <- Inf
    min.i <- which.min(dist.i)
    point.new <- labels[min.i, 1:2]
    labels[min.i, "flag"] <- TRUE
    
    if (count > 0) {
      mesh.num <- as.integer(dist.i[min.i])
      mesh.step <- seq(1e-2, 1 - 1e-2, 1 / mesh.num)
      mesh.i <- data.frame(
        x = (point.old$x * mesh.step) + (point.new$x * (1 - mesh.step)),
        y = (point.old$y * mesh.step) + (point.new$y * (1 - mesh.step))
      )
      labels.mesh <- c(labels.mesh, list(mesh.i))
    }
    labels.mesh <- c(labels.mesh, list(point.new))
    point.old <- point.new
    count <- count + 1
  }
  
  labels.mesh.df <- do.call("rbind", labels.mesh)
  labels.mesh.df
}

create.all.mesh <- function(label.df) {
  label.PE <- create.label.mesh(
    labels = label.df[label.df$label == "PE", 2:3],
    start = label.df[label.df$label == "START", 2:3]
  )
  label.B1 <- create.label.mesh(
    labels = label.df[label.df$label == "B1", 2:3],
    start = label.df[label.df$label == "START", 2:3]
  )
  label.B2 <- create.label.mesh(
    labels = label.df[label.df$label == "B2", 2:3],
    start = label.df[label.df$label == "START", 2:3]
  )
  
  list(PE = label.PE, B1 = label.B1, B2 = label.B2)
}

# Create empty retina background ####
create.empty.bk <- function(pic, x0, y0) {
  pic.bk <- pic +
    geom_circle(
      data = data.frame(),
      aes(x0 = x0, y0 = y0, r = 7),
      fill = "#BBBBBB", color = "#777777"
    ) +
    geom_circle(
      aes(x0 = x0, y0 = y0, r = 6),
      fill = "#DDDDDD", color = "#777777"
    ) +
    geom_circle(
      aes(x0 = x0, y0 = y0, r = 5),
      fill = "#FFFFFF", color = "#777777"
    ) +
    geom_ellipse(
      aes(x0 = x0 + 5.2, y0 = y0, a = 2.4, b = 4.0, angle = 0),
      fill = "#FFFFFF", color = "#777777"
    ) + new_scale_color()
  
  pic.bk
}

# DEG number at different resolutions ####
diff.res.num <- function(obj, min.res, max.res, res.gap) {
  deg.num.pos <- c()
  deg.num.neg <- c()
  res.seq <- seq(min.res, max.res, by = res.gap)
  
  for (res.i in res.seq) {
    obj <- SetIdent(obj, value = paste0("Spatial_mix_res.", res.i))
    obj.diff <- FindAllMarkers(
      obj, test.use = "wilcox",
      logfc.threshold = 0.25, min.pct = 0.1
    )
    obj.diff <- obj.diff[obj.diff$p_val_adj < 0.05, ]
    
    # Calculate number
    obj.diff.pos <- obj.diff %>%
      filter(avg_log2FC > 0) %>%
      count(cluster, name = "number")
    if (nrow(obj.diff.pos) == nlevels(obj)) {
      obj.diff.pos.num <- min(obj.diff.pos$number)
    } else {
      obj.diff.pos.num <- 0
    }
    deg.num.pos[as.character(res.i)] <- obj.diff.pos.num
    
    obj.diff.neg <- obj.diff %>%
      filter(avg_log2FC < 0) %>%
      count(cluster, name = "number")
    if (nrow(obj.diff.neg) == nlevels(obj)) {
      obj.diff.neg.num <- min(obj.diff.neg$number)
    } else {
      obj.diff.neg.num <- 0
    }
    deg.num.neg[as.character(res.i)] <- obj.diff.neg.num
  }
  
  list(pos = deg.num.pos, neg = deg.num.neg)
}
