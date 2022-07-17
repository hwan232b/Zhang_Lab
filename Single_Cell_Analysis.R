suppressPackageStartupMessages({
  # imports for analyses
  library(symphony)
  library(Seurat)
  library(dplyr)
  library(singlecellmethods)
  library(harmony)
  library(irlba)
  library(gridExtra)
  library(tidyverse)
  library(msigdbr)
  library(parallel)
  library(magrittr)
  
  # imports for figures
  library(viridis)
  library(ggrepel)
  library(ggrastr)
  library(ggpubr)
  library(ggplot2)
  library(ggthemes)
  library(ggpointdensity)
  library(cowplot)
  
  # linear modeling
  library(nlme)
  library(limma)
  library(glmnet)
  library(stringr) 
})

fig.size <- function (height, width) {
  options(repr.plot.height = height, repr.plot.width = width)
}

meta_colors <- list(
  "cluster" = c(
    "SC-F1" = "#6BAED6",
    "SC-F2" = "#08306B", 
    "SC-F3" = "#DEEBF7",
    "SC-F4" = "grey",
    "SC-T1" = "#8C510A",
    "SC-T2" = "brown",
    "SC-T3" = "#FFFF33",
    "SC-T4" = "#C7EAE5",
    "SC-T5" = "#003C30",
    "SC-T6" = "#35978F", 
    "SC-B1" = "#FCBBA1",
    "SC-B2" = "#CB181D", 
    "SC-B3" = "#67000D",
    "SC-B4" = "#FB9A99",
    "SC-M1" = "#AE017E",
    "SC-M2" = "#F768A1",
    "SC-M3" = "#FDE0EF", 
    "SC-M4" = "#49006A",
    "Centroid" = "grey"
  )
)

# Functions
FindVariableGenesBatch <- function(exprs_mat, meta_df, genes_exclude = NULL, ngenes_use = 1e3, expr_min = .1) {
  if (!is.null(genes_exclude)) {
    genes_use <- setdiff(row.names(exprs_mat), genes_exclude)
  }
  else #hannah bug fix
  {
    genes_use <- row.names(exprs_mat)
  }
  x_res <- split(meta_df$cell, meta_df$sample) %>% lapply(function(x) {
    FindVariableGenesSeurat(exprs_mat[genes_use, x]) %>% 
      subset(gene.mean >= expr_min) %>% 
      tibble::rownames_to_column("gene") %>% 
      dplyr::arrange(-gene.dispersion) %>%
      head(ngenes_use)
  })
  data.table(Reduce(rbind, x_res))[, .N, by = gene][order(-N)]    
}


FindVariableGenesSeurat <- function (data, x.low.cutoff = 0.1, x.high.cutoff = 8,
                                     y.cutoff = 1, y.high.cutoff = Inf, num.bin = 0,
                                     binning.method = "equal_width", sort.results = TRUE,
                                     display.progress = TRUE, ...)
{
  genes.use <- rownames(data)
  if (class(data) != "dgCMatrix") {
    data <- as(as.matrix(data), "dgCMatrix")
  }
  ## (1) get means and variances
  gene.mean <- FastExpMean(data, display.progress)
  names(gene.mean) <- genes.use
  gene.dispersion <- FastLogVMR(data, display.progress)
  names(gene.dispersion) <- genes.use
  
  gene.dispersion[is.na(x = gene.dispersion)] <- 0
  gene.mean[is.na(x = gene.mean)] <- 0
  mv.df <- data.frame(gene.mean, gene.dispersion)
  rownames(mv.df) <- rownames(data)
  
  ## (OPTIONAL) do the binning correction
  if (num.bin > 0) {
    if (binning.method == "equal_width") {
      data_x_bin <- cut(x = gene.mean, breaks = num.bin)
    }
    else if (binning.method == "equal_frequency") {
      data_x_bin <- cut(x = gene.mean, breaks = c(-1, quantile(gene.mean[gene.mean >
                                                                           0], probs = seq(0, 1, length.out = num.bin))))
    }
    else {
      stop(paste0("Invalid selection: '", binning.method,
                  "' for 'binning.method'."))
    }
    names(x = data_x_bin) <- names(x = gene.mean)
    mean_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                     FUN = mean)
    sd_y <- tapply(X = gene.dispersion, INDEX = data_x_bin,
                   FUN = sd)
    gene.dispersion.scaled <- (gene.dispersion - mean_y[as.numeric(x = data_x_bin)])/sd_y[as.numeric(x = data_x_bin)]
    gene.dispersion.scaled[is.na(x = gene.dispersion.scaled)] <- 0
    ##names(gene.dispersion.scaled) <- names(gene.mean)
    
    mv.df$gene.dispersion.scaled <- gene.dispersion.scaled
  }
  
  return(mv.df)
}

environment(FindVariableGenesSeurat) <- asNamespace("Seurat")

ScaleDataSeurat <- function (data.use, margin = 1, scale.max = 10,
                             block.size = 1000) {
  
  if (margin == 2) data.use %<>% t
  max.block <- ceiling(nrow(data.use)/block.size)
  
  ## Define data and functions to use in sparse and dense cases
  if (class(data.use) == "dgCMatrix" | class(data.use) == "dgTMatrix") {
    scale_fxn <- function(x) {
      FastSparseRowScale(mat = x, scale = TRUE, center = TRUE,
                         scale_max = scale.max, display_progress = FALSE)
    }
  } else {
    scale_fxn <- function(x) {
      FastRowScale(mat = x, scale = TRUE, center = TRUE,
                   scale_max = scale.max, display_progress = FALSE)
    }
    data.use <- as.matrix(data.use)
  }
  
  ## Do scaling, at once or in chunks
  if (max.block == 1) {
    scaled.data <- scale_fxn(data.use)
  } else {
    scaled.data <- matrix(NA, nrow(data.use), ncol(data.use))
    for (i in 1:max.block) {
      idx.min <- (block.size * (i - 1))
      idx.max <- min(nrow(data.use), (block.size * i - 1) + 1)
      my.inds <- idx.min:idx.max
      scaled.data[my.inds, ] <- scale_fxn(data.use[my.inds, , drop = F])
    }
  }
  
  colnames(scaled.data) <- colnames(data.use)
  row.names(scaled.data) <- row.names(data.use)
  scaled.data[is.na(scaled.data)] <- 0
  if (margin == 2) scaled.data %<>% t
  return(scaled.data)
}
environment(ScaleDataSeurat) <- asNamespace("Seurat")


fig.size <- function(height, width) {
  options(repr.plot.height = height, repr.plot.width = width)
}

SingleFeaturePlotSeurat <- function (data.use, feature, data.plot, pt.size, pch.use, cols.use,
                                     dim.codes, min.cutoff, max.cutoff, coord.fixed, no.axes,
                                     no.title = FALSE, no.legend, dark.theme, vector.friendly = FALSE,
                                     png.file = NULL, png.arguments = c(10, 10, 100))
{
  if (vector.friendly) {
    previous_call <- blank_call <- png_call <- match.call()
    blank_call$pt.size <- -1
    blank_call$vector.friendly <- FALSE
    png_call$no.axes <- TRUE
    png_call$no.legend <- TRUE
    png_call$vector.friendly <- FALSE
    png_call$no.title <- TRUE
    blank_plot <- eval(blank_call, sys.frame(sys.parent()))
    png_plot <- eval(png_call, sys.frame(sys.parent()))
    png.file <- SetIfNull(x = png.file, default = paste0(tempfile(),
                                                         ".png"))
    ggsave(filename = png.file, plot = png_plot, width = png.arguments[1],
           height = png.arguments[2], dpi = png.arguments[3])
    to_return <- AugmentPlot(blank_plot, png.file)
    file.remove(png.file)
    return(to_return)
  }
  idx.keep <- which(!is.na(data.use[feature, ]))
  data.gene <- data.frame(data.use[feature, idx.keep])
  #     data.gene <- na.omit(object = data.frame(data.use[feature,
  #         ]))
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.gene)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.gene)
  data.gene <- sapply(X = data.gene, FUN = function(x) {
    return(ifelse(test = x < min.cutoff, yes = min.cutoff,
                  no = x))
  })
  data.gene <- sapply(X = data.gene, FUN = function(x) {
    return(ifelse(test = x > max.cutoff, yes = max.cutoff,
                  no = x))
  })
  data_plot <- data.plot[idx.keep, ]
  data_plot$gene <- data.gene
  if (length(x = cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  }
  else {
    brewer.gran <- length(x = cols.use)
  }
  if (all(data.gene == 0)) {
    data.cut <- 0
  }
  else {
    data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene),
                                                 breaks = brewer.gran)))
  }
  data_plot$col <- as.factor(x = data.cut)
  p <- data_plot %>%
    dplyr::arrange(col) %>%
    ggplot(mapping = aes(x = x, y = y))
  if (brewer.gran != 2) {
    if (length(x = cols.use) == 1) {
      p <- p + geom_point(mapping = aes(color = col), size = pt.size,
                          shape = pch.use) + #scale_color_brewer(palette = cols.use)
        scale_color_viridis(option = "plasma", end = .9)
    }
    else {
      p <- p + geom_point(mapping = aes(color = col), size = pt.size,
                          shape = pch.use) + #scale_color_manual(values = cols.use)
        scale_color_viridis(option = "plasma", end = .9)
    }
  }
  else {
    if (all(data_plot$gene == data_plot$gene[1])) {
      warning(paste0("All cells have the same value of ",
                     feature, "."))
      p <- p + geom_point(color = cols.use[1], size = pt.size,
                          shape = pch.use)
    }
    else {
      p <- p + geom_point(mapping = aes(color = gene),
                          size = pt.size, shape = pch.use) + scale_color_viridis(option = "plasma", end = .9
                          )
    }
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  if (no.axes) {
    p <- p + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                   axis.text.y = element_blank(), axis.ticks = element_blank(),
                   axis.title.x = element_blank(), axis.title.y = element_blank())
    if (!no.title)
      p <- p + labs(title = feature, x = "", y = "")
    if (no.title)
      p <- p + labs(x = "", y = "")
  }
  else {
    if (no.title)
      p <- p + labs(x = dim.codes[1], y = dim.codes[2])
    if (!(no.title))
      p <- p + labs(title = feature) + labs(x = "", y = "")
  }
  if (no.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (coord.fixed) {
    p <- p + coord_fixed()
  }
  return(p)
}
environment(SingleFeaturePlotSeurat) <- asNamespace("Seurat")

PlotFeatures <- function(umap_use, features_plot, exprs_use, cells_use, ncols, pt_size = .5, pt_shape = ".", q_lo = "q10", q_hi = "q90") {
  if (missing(cells_use)) cells_use <- 1:nrow(umap_use)
  if (missing(ncols)) ncols <- round(sqrt(length(features_plot)))
  
  plt_list <- lapply(features_plot, function(feature_use) {
    SingleFeaturePlotSeurat(exprs_use[, cells_use], feature_use, data.frame(x = umap_use[cells_use, 1], y = umap_use[cells_use, 2]),
                            pt.size = pt_size, pch.use = pt_shape, cols.use = c("lightgrey", "blue"),
                            dim.codes = c("UMAP 1", "UMAP 2"), min.cutoff = c(q10 = q_lo), max.cutoff = c(q90 = q_hi),
                            coord.fixed = FALSE, no.axes = FALSE, dark.theme = FALSE, no.legend = TRUE)
  })
  plot_grid(plotlist = plt_list, ncol = ncols)
  #return(plt_list)
}

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0) {
  my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
  nn.ranked <- my.knn$nn.idx
  
  snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
  rownames(snn_res) <- row.names(data.use)
  colnames(snn_res) <- row.names(data.use)
  return(snn_res)
}
environment(BuildSNNSeurat) <- asNamespace("Seurat")

NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- scaling_factor * A@x
  if (do_ftt) {
    A@x <- sqrt(A@x) + sqrt(1 + A@x)
  } else {
    A@x <- log(1 + A@x)
  }
  return(A)
}

plot_clusters3 <- function(cluster_ids, labels, pt_size = 14, umap_use = umap_post, do_labels = FALSE) {
  cluster_table <- table(cluster_ids)
  clusters_keep <- names(which(cluster_table > 20))
  plt_df <- umap_use %>% data.frame() %>% cbind(cluster = cluster_ids) %>%
    subset(cluster %in% clusters_keep) 
  plt <- plt_df %>% 
    ggplot(aes(X1, X2, col = factor(cluster))) + geom_point(shape = '.', alpha = .6) + 
    theme_tufte() + geom_rangeframe(col = "black") + 
    #         theme(axis.line = element_line()) +
    guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 21, size = 4))) + 
    scale_color_manual(values = singler.colors) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(hjust = .5)) + 
    guides(col = FALSE)
  
  if (do_labels) 
    plt <- plt + geom_label(data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = cluster], 
                            aes(label = cluster), size = pt_size, alpha = .8)
  return(plt)
}

# Load AMP phase 1 RA single-cell data

crohns_exprs_norm <- readRDS("/Users/hannahwang/github/Projects/Zhang_lab/crohns_exprs_norm_qc_kuhn.rds")
crohns_meta_all <- readRDS("/Users/hannahwang/github/Projects/Zhang_lab/crohns_meta_kuhn.rds")

HIV_exprs_norm <- readRDS("/Users/hannahwang/github/Projects/Zhang_lab/HIV_exprs_norm_qc_ShaoboWang.rds")
HIV_meta_all <- readRDS("/Users/hannahwang/github/Projects/Zhang_lab/HIV_meta_ShaoboWang.rds")

length(intersect(row.names(crohns_exprs_norm),row.names(HIV_exprs_norm)))

intersect1 <- intersect(row.names(crohns_exprs_norm),row.names(HIV_exprs_norm))

HIV_intersect <- HIV_exprs_norm[intersect1,]
crohns_intersect <- crohns_exprs_norm[intersect1,]

# add cluster data into crohns_meta_all
# crohns_combined_meta <- read.csv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_All.combined.metadata.csv.gz")
# crohns_combined_meta1 <- sub("_", "-", crohns_combined_meta[,1])
# rownames(crohns_combined_meta) <- crohns_combined_meta1
# crohns_combined_meta1<-as.data.frame(crohns_combined_meta1)
# 
# crohns_meta <- readRDS("crohns_kuhn.rds")
# rownames(crohns_meta) <- sub(".*_ ", "", rownames(crohns_meta))
# crohns_meta <- cbind(crohns_meta,rownames(crohns_meta))
# matrix_crohns <- matrix(nrow=62166, ncol=6)
# df_crohns <- as.data.frame(matrix_crohns)
# crohns_meta <- cbind(crohns_meta,df_crohns)
# 
# new_meta <- matrix(nrow=62166, ncol=2)
# 
# for (row in 1:nrow(crohns_meta)) {
#   if((crohns_combined_meta1[row,1] %in% crohns_meta[,4])){
#     new_meta[row,2] <- (crohns_meta[row,1])
#     new_meta[row,2] <- (crohns_combined_meta[row,9])
#   } else {crohns_meta[row,4] <- (NA)
#   }} 
# 
# crohns_meta_all <- cbind(crohns_meta_all,crohns_meta[,9])
# colnames(crohns_meta_all)[7] <- "cluster"

# new code
crohns_combined_meta <- read.csv("/Users/hannahwang/github/Projects/Zhang_lab/GSE163314_All.combined.metadata.csv.gz")
crohns_combined_meta[,2] <- sub("_.*", "",crohns_combined_meta[,2])
crohns_combined_meta[,2] <- paste0(crohns_combined_meta[,2],"_")
crohns_combined_meta[,1] <- paste(crohns_combined_meta[,2],crohns_combined_meta[,1])
crohns_combined_meta[,1] <- sub("_", "-", crohns_combined_meta[,1])
crohns_combined_meta[,1] <- sub("_", "-", crohns_combined_meta[,1])
crohns_combined_meta[,1] <- sub("- ", "_", crohns_combined_meta[,1])
crohns_combined_meta[,1] <- sub("-.*", "",crohns_combined_meta[,1])


#crohns_meta <- readRDS("crohns_kuhn.rds")
# X_df <- as.data.frame(matrix(rownames(crohns_meta)))
# crohns_meta <- cbind(crohns_meta,X_df)
# crohns_meta[,4] <- sub("P", "B", crohns_meta[,4])
# crohns_meta[,4] <- sub("_", "-", crohns_meta[,4])
# crohns_meta[,4] <- sub("- ", "_", crohns_meta[,4])
# crohns_meta[,4] <- sub("-.*", "", crohns_meta[,4])

# 25,000 blood samples were mislabeled with P
crohns_meta_all[,1] <- sub("P", "B", crohns_meta_all[,1])

# our query to verify that these are in fact blood tissues
p1_df = crohns_meta_all[grepl('^P', crohns_meta_all$sample) & crohns_meta_all$tissue=="Blood",]
p2_df = crohns_meta_all[grepl('^P', crohns_meta_all$sample) & crohns_meta_all$tissue!="Blood",]

crohns_meta_all[,1] <- sub("-.*", "", crohns_meta_all[,1])
crohns_meta_all[,1] <- sub(" ", "", crohns_meta_all[,1])

#merge_combined <- merge(x=crohns_meta,y=crohns_combined_meta,by.x="V1", by.y="X",all.x=TRUE)
# clusters1 <- cbind(merge_combined$V1,merge_combined$clusters)
# clusters1<- as.data.frame(clusters1)
# colnames(clusters1) <- c("V1","clusters")
#clusters1 <- merge_combined[c("V1","clusters")]

crohns_meta_all <- merge(x=crohns_meta_all,y=crohns_combined_meta,by.x="cell", by.y="X",all.x=TRUE)
crohns_meta_all <- crohns_meta_all[c("cell", "sample", "disease", "batch", "tissue", "dataset", "clusters")]

# add empty cluster column to HIV meta
HIV_empty <- matrix(nrow=35750,ncol=1)
colnames(HIV_empty)[1] <- "clusters"
HIV_empty <- as.data.frame(HIV_empty)
HIV_meta_all <- cbind(HIV_meta_all,HIV_empty)
HIV_meta_all[,1] <- sub("-.*", "", HIV_meta_all[,1])
HIV_meta_all[,1] <- sub(" ", "", HIV_meta_all[,1])

# crohns_intersect <- crohns_exprs_norm[rownames(crohns_exprs_norm) %in% intersect1, ]
# HIV_intersect <- HIV_exprs_norm[rownames(HIV_exprs_norm) %in% intersect1, ]
# NA <- c("NA")
# crohns_intersect <- crohns_intersect[!rownames(crohns_intersect) %in% NA, ]
# HIV_intersect <- HIV_intersect[!rownames(HIV_intersect) %in% NA, ]
# crohns_intersect <- unique(rownames(crohns_intersect))
# HIV_intersect <- unique(rownames(HIV_intersect))

exprs_norm <- cbind(crohns_intersect,HIV_intersect)
meta_all <- rbind(HIV_meta_all,crohns_meta_all)

#intersection of meta_all$cell and exprs_norm to make them exact matches
colnames(exprs_norm) <- sub("-.*", "", colnames(exprs_norm))
colnames(exprs_norm) <- sub(" ", "", colnames(exprs_norm))


intersect2 <- intersect(colnames(exprs_norm),meta_all$cell)
meta_all <- meta_all[meta_all$cell %in% intersect2,]

# function "all" order matters so we are alphabetizing
exprs_norm <- exprs_norm[,order(colnames(exprs_norm))]
meta_all <- meta_all[order(rownames(meta_all)),]

# check if TRUE
all(colnames(exprs_norm) == meta_all$cell)

# removing all rownames with missing values 
exprs_norm <- exprs_norm[-2,]

# finding highly variable sample
genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-", row.names(exprs_norm), value = TRUE)
vargenes_df <- FindVariableGenesBatch(exprs_norm, meta_all)

nrow(vargenes_df)
var_genes <- vargenes_df$gene

exprs_norm <- exprs_norm[var_genes,]

exprs_norm <- Matrix::Matrix(exprs_norm, sparse = TRUE)
class(exprs_norm)

vargenes_means_sds <- tibble(symbol = var_genes, mean = Matrix::rowMeans(exprs_norm))
vargenes_means_sds$stddev <- singlecellmethods::rowSDs(exprs_norm, vargenes_means_sds$mean)

vargenes_means_sds[1:4,]
dim(vargenes_means_sds)

# make graph
options(repr.plot.height = 5, repr.plot.width = 7)
ggplot(vargenes_means_sds, aes(mean, stddev)) +
  geom_pointdensity(size = 1) +
  scale_color_viridis() +
  theme_bw(base_size = 20) +
  theme(legend.position="none")

# Scale data
ref_exp_scaled <- singlecellmethods::scaleDataWithStats(exprs_norm, vargenes_means_sds$mean, 
                                                        vargenes_means_sds$stddev, 1)

# Run SVD, save gene loadings
s = irlba(ref_exp_scaled, nv = 20)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

# Gene loadings: genes x 20
dim(loadings)
loadings[1:4,]

# cell loadings: 20 x cell number
Z_pca_ref[1:4, ]

# Run Harmony
ref_harmObj = harmony::HarmonyMatrix(
  data_mat = t(Z_pca_ref), ## PCA embedding matrix of cells
  meta_data = meta_all, ## dataframe with cell labels
  theta = c(2), ## cluster diversity enforcement
  vars_use = c('sample'), ## variable to integrate out
  nclust = 200, ## number of clusters in Harmony model: use more cluster centroids to capture more subsets
  max.iter.harmony = 10,
  return_object = TRUE, ## return the full Harmony model object
  do_pca = FALSE ## don't recompute PCs
)

# buildReferenceFromHarmonyObj()
# Compress Harmony reference into a Symphony reference
reference = symphony::buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  meta_all,
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs
  verbose = TRUE,
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  )

# save
saveRDS(reference, '2022_07_14_single_cell_analysis_clusters.rds')


# Visualize
umap_labels <- cbind(reference$umap, reference$meta_data)
colnames(umap_labels)[1:2] <- c("UMAP1", "UMAP2")
umap_labels[1:4,]

str(reference)

# define the colors for each plot
library(RColorBrewer)

meta_colors <- list(
  "disease" = c(
    "AS" = "#EFFF03",
    "CD" = "green", 
    "CDAS" = "purple",
    "Control" = "blue",
    "Healthy" = "orange",
    "HIV" = "red"
  ))

meta_colors <- list(
  "clusters" = c(
    "Bmem" = "purple",
    "Bnaive" = "purple", 
    "CD4_Tmem" = "purple",
    "CD4_Tnaive" = "purple",
    "CD8_T" = "purple",
    "DC2_CD1C" = "purple",
    "DC4_CD16"= "purple",
    "Mono1_CD14" = "purple",
    "NK" = "purple",
    "Plasma" = "purple",
    "SG2M" = "purple",
    "Tcell" = "purple",
    "DC6_pDC" = "purple"
  ))

meta_colors <- list(
  "tissue" = c(
    "blood" = "orange",
    "Blood" = "green", 
    "Colon" = "purple"
  ))

p1 <- ggplot(umap_labels[sample(nrow(umap_labels)),],
             aes(x = UMAP1, y = UMAP2, fill= sample)
) +
 # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  #facet_wrap(~sample)+
  #scale_fill_manual(values = meta_colors$cluster, name = "") +
  theme_bw(base_size = 15) 


p2 <- ggplot(umap_labels[sample(nrow(umap_labels)),],
             aes(x = UMAP1, y = UMAP2, fill= disease)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  #facet_grid(disease ~ tissue) +
  #facet_wrap(~disease)+
  scale_fill_manual(values = meta_colors$disease, name = "") +
  theme_bw(base_size = 15) 

p3 <- ggplot(umap_labels[sample(nrow(umap_labels)),],
             aes(x = UMAP1, y = UMAP2, fill= batch)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  # scale_fill_manual(values = meta_colors$cluster, name = "") +
  theme_bw(base_size = 15) 

p4 <- ggplot(umap_labels[sample(nrow(umap_labels)),],
             aes(x = UMAP1, y = UMAP2, fill= tissue)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  scale_fill_manual(values = meta_colors$tissue, name = "") +
  #facet_wrap(~tissue)+
  theme_bw(base_size = 15) 

p5 <- ggplot(umap_labels[sample(nrow(umap_labels)),],
             aes(x = UMAP1, y = UMAP2, fill= dataset)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  # scale_fill_manual(values = meta_colors$cluster, name = "") +
  #facet_wrap(~dataset)
  theme_bw(base_size = 15) 

umap_labels1 <- umap_labels[complete.cases(umap_labels), ]      

p6 <- ggplot(umap_labels1,
             aes(x = UMAP1, y = UMAP2, fill= clusters)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  #scale_fill_manual(values = meta_colors$clusters, name = "") +
  theme_bw(base_size = 15) +
  facet_wrap(~clusters, nrow =4)
  #scale_fill_brewer(palette="Set3")

meta_colors <- list(
  "clusters1" = c(
    "Bmem" = "#e41a1c",
    "Bnaive" = "#377eb8", 
    "CD4_Tmem" = "#4daf4a",
    "CD4_Tnaive" = "purple",
    "CD8_T" = "orange",
    "DC2_CD1C" = "#ec34be",
    "DC4_CD16"= "#09c6ce",
    "Mono1_CD14" = "#633126",
    "NK" = "#167140",
    "Plasma" = "black",
    "SG2M" = "#c7db1a",
    "Tcell" = "#7d7e69",
    "DC6_pDC" = "#1c3ea5"
  ))

p7 <- ggplot(umap_labels1,
             aes(x = UMAP1, y = UMAP2, fill= clusters)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  #scale_fill_manual(values = meta_colors$clusters1, name = "") +
  theme_bw(base_size = 15)
  #facet_wrap(~clusters, nrow=4)
#scale_fill_brewer(palette="Set3")

options(repr.plot.height = 5, repr.plot.width = 15)
plot_grid(p1, p2, p4, p5, labels = c('A', 'B', 'C', 'D'), ncol =2)
plot_grid(p6,p7)

dim(Z_pca_ref)
dim(umap_labels)
pcs <- as.data.frame(t(Z_pca_ref))
colnames(pcs) <- paste0("PC", colnames(pcs), sep="")
umap_labels_pcs <- cbind(umap_labels, pcs)
umap_labels_pcs[1:4,]
x1 <- ggplot(umap_labels_pcs[sample(nrow(umap_labels_pcs)),],
             aes(x = PCV1, y = PCV3, fill= sample)
) +
  # geom_hex(bins = 150) +
  geom_point(size = 1, stroke = 0.0001, shape = 21, alpha = 0.6) +
  #scale_fill_manual(values = meta_colors$cluster, name = "") +
  theme_bw(base_size = 15) 
x1



# Get centroid
cluster_sizes = reference$cache[[1]] %>% as.matrix()
centroid_sums = t(reference$Z_corr %*% t(reference$R)) %>% as.data.frame()
centroids = sweep(centroid_sums, 1, cluster_sizes, "/")
colnames(centroids) = paste0("hPC", c(1:20))
dim(centroids)

ref_umap_model = uwot::load_uwot(reference$save_uwot_path, verbose = FALSE)
umap_centroids = uwot::umap_transform(centroids, ref_umap_model)
umap_centroids <- as.data.frame(umap_centroids)
colnames(umap_centroids) <- c("UMAP1", "UMAP2")
umap_centroids$cell <- rep("NA", nrow(umap_centroids))
umap_centroids$cell_type <- rep("NA", nrow(umap_centroids))
umap_centroids$disease <- rep("NA", nrow(umap_centroids))
umap_centroids$sample <- rep("NA", nrow(umap_centroids))
umap_centroids$plate <- rep("NA", nrow(umap_centroids))
umap_centroids$cluster <- rep("Centroid", nrow(umap_centroids))
umap_centroids[1:4,]
dim(umap_centroids)
