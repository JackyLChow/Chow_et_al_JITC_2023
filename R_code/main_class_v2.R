# load packages
library(SingleCellExperiment)
library(scater)

sce <- readRDS("~/_Projects/Chow Nat Comm/data/sce_FZ02.RDS")
sce <- logNormCounts(sce)

sce$imm_call <- NA
sce$class <- NA

set.seed(415); sce <- runUMAP(sce)

sce$imm_call <- ifelse(reducedDim(sce, "UMAP")[, 2] > -5, "CD45hi", "CD45lo")

sce$class <- ifelse(sce$imm_call == "CD45hi" & reducedDim(sce, "UMAP")[, 1] > -5, "tnk",
                    ifelse(sce$imm_call == "CD45hi" & reducedDim(sce, "UMAP")[, 1] < -5, "mye", NA))

nim <- sce[, sce$imm_call == "CD45lo"]

set.seed(415); nim <- runUMAP(nim)

nim$class <- ifelse(reducedDim(nim, "UMAP")[, 1] < -10, "mst",
                    ifelse(reducedDim(nim, "UMAP")[, 2] > 0, "tum",
                           ifelse(reducedDim(nim, "UMAP")[, 2] > -6, "ves", "b_c")))

sce[, colnames(nim)]$class <- nim$class

#saveRDS(sce, "~/_Projects/Chow Nat Comm/data/sce_classed_FZ02.RDS")
#sce <- readRDS("~/_Projects/Chow Nat Comm/data/sce_classed_FZ02.RDS")

tnk <- sce[, sce$class == "tnk"]

set.seed(415); tnk <- runUMAP(tnk)
set.seed(415); tnk <- runTSNE(tnk)

hc <- hclust(dist(reducedDim(tnk, "UMAP")))
tnk$hc_14 <- factor(cutree(hc, 14))

#saveRDS(tnk, "~/_Projects/Chow Nat Comm/data/tnk_exp.RDS")

mye <- sce[, sce$class == "mye"]

set.seed(415); mye <- runUMAP(mye)
set.seed(415); mye <- runTSNE(mye)

hc <- hclust(dist(reducedDim(mye, "UMAP")))
mye$hc_16 <- factor(cutree(hc, 16))

#saveRDS(mye, "~/_Projects/Chow Nat Comm/data/mye_exp.RDS")

for(i in c("b_c", "tum", "mst", "ves")){
  sce_ <- sce[, sce$class == i]
  saveRDS(sce_, paste0("~/_Projects/Chow Nat Comm/data/", i, "_FZ02.RDS"))
}















