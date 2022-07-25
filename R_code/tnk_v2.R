# load packages
library(gridExtra)
library(ggplotify)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(circlize)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(readr)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(tradeSeq)
library(limma)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

# load data
tnk <- readRDS("~/_Projects/Chow Nat Comm/data/tnk_exp.RDS")
# load human gene labeling database
hs <- org.Hs.eg.db

##########~~~TNK cluster subclassification~~~##########
# generate plot data.frame using TSNE coordinates
pt <- cbind(x = reducedDim(tnk, "TSNE")[, 1],
            y = reducedDim(tnk, "TSNE")[, 2],
            data.frame(t(logcounts(tnk)[c("CD3D", "CD8A", "CD4", "KLRD1",
                                          "TCF7", "IL7R", "HAVCR2", "TOX",
                                          "FOXP3", "IL2RA", "IFNG", "MKI67"), ])),
            data.frame(colData(tnk)))

# cell subclassification
## assign subclass after examining variant and hallmark gene expression
pt$subclass <- ifelse(pt$hc_14 == 1, "cd4_3",
             ifelse(pt$hc_14 == 2, "cd4_1",
             ifelse(pt$hc_14 == 3, "cd8_1",
             ifelse(pt$hc_14 == 4, "cd8_5",
             ifelse(pt$hc_14 == 5, "n_k_1",
             ifelse(pt$hc_14 == 6, "cd4_2",
             ifelse(pt$hc_14 == 7, "cd4_4",
             ifelse(pt$hc_14 == 8, "n_k_2",
             ifelse(pt$hc_14 == 9, "cd8_4",
             ifelse(pt$hc_14 == 10, "cd8_2",
             ifelse(pt$hc_14 == 11, "cd8_3",
             ifelse(pt$hc_14 == 12, "nkt_3",
             ifelse(pt$hc_14 == 13, "nkt_1",
             ifelse(pt$hc_14 == 14, "nkt_2", NA
             ))))))))))))))

sub_nom <- data.frame(subclass = c(paste0("cd4_", 1:4),
                                   paste0("cd8_", 1:5),
                                   paste0("n_k_", 1:2),
                                   paste0("nkt_", 1:3)),
                      subclass_fine = c("CD4_MALAT1", "CD4_FOS", "CD4_ANXA1", "CD4_FOXP3",
                                        "CD8_TCF7", "CD8_ITM2C", "CD8_CCR5", "CD8_HAVCR2", "CD8_MKI67",
                                        "NK_NCAM1", "NK_FCGR3A",
                                        "NKT_KLRG1", "NKT_GZMH", "NKT_GZMK"))

pt <- left_join(pt, sub_nom, by = "subclass")

tnk$subclass <- pt$subclass
tnk$subclass_fine <- pt$subclass_fine

# saveRDS(tnk, "~/_Projects/Chow Nat Comm/data/tnk_FZ02.RDS")
# tnk <- readRDS("~/_Projects/Chow Nat Comm/data/tnk_FZ02.RDS")

sce <- SingleCellExperiment()
for(file_ in list.files("~/_Projects/Chow Nat Comm/data/sce_final", full.names = T)){
  if(nrow(sce) == 0){
    sce <- readRDS(file_)
  } else {
    sce <- cbind(sce, readRDS(file_))
  }
}

tnk <- sce[, sce$class == "tnk"]
reducedDim(tnk, "UMAP") <- readRDS("~/_Projects/Chow Nat Comm/data/tnk_UMAP.RDS")

# generate plot data.frame using TSNE coordinates
pt <- cbind(x = reducedDim(tnk, "TSNE")[, 1],
            y = reducedDim(tnk, "TSNE")[, 2],
            data.frame(t(logcounts(tnk)[c("CD3D", "CD8A", "CD4", "KLRD1",
                                          "TCF7", "IL7R", "HAVCR2", "TOX",
                                          "FOXP3", "IL2RA", "IFNG", "MKI67"), ])),
            data.frame(colData(tnk)))

sc_col <- colorRampPalette(c("#8d5103",
                             "#f0f0c0",
                             "#f76f05",
                             "#ffd600"))(14)

## plot on parent TSNE
########## by subclass ##########
labels <- data.frame(label = c("CD4_MALAT1", "CD4_FOS", "CD4_ANXA1", "CD4_FOXP3",
                               "CD8_TCF7", "CD8_ITM2C", "CD8_CCR5", "CD8_HAVCR2", "CD8_MKI67",
                               "NK_NCAM1", "NK_FCGR3A",
                               "NKT_KLRG1", "NKT_GZMH", "NKT_GZMK"),
                     x = c(-12, 8, 20, -5,
                           2, -12, -15, -28, -29,
                           21, 12,
                           0, -9, -2),
                     y = c(-6, -6, 5, -15,
                           12, 16, -2, 15, -5,
                           23, 33,
                           32, 27, 15))

png(filename = "Fig2_TNSE_sc.png", width = 1250, height = 1250, res = 300)

ggplot(pt, aes(x, y)) +
  geom_point(size = 0.1,
             aes(color = subclass)) +
  geom_text_repel(data = filter(labels,
                                grepl("KLRG", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 10,
                  nudge_y = 12) +
  geom_text_repel(data = filter(labels,
                                grepl("GZMH", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 4,
                  nudge_y = 12) +
  geom_text_repel(data = filter(labels,
                                grepl("FCGR|NCAM", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 10,
                  nudge_y = 7) +
  geom_text_repel(data = filter(labels,
                                grepl("TCF7|GZMK", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 33,
                  nudge_y = 3) +
  geom_text_repel(data = filter(labels,
                                grepl("ANXA", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 10,
                  nudge_y = 3) +
  geom_text_repel(data = filter(labels,
                                grepl("FOS|FOXP3", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 5,
                  nudge_y = -5) +
  geom_text_repel(data = filter(labels,
                                grepl("CCR5|MALA", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -10,
                  nudge_y = -15) +
  geom_text_repel(data = filter(labels,
                                grepl("MKI", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -20,
                  nudge_y = -5) +
  geom_text_repel(data = filter(labels,
                                grepl("HAVC", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -10,
                  nudge_y = 6) +
  geom_text_repel(data = filter(labels,
                                grepl("ITM2", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -10,
                  nudge_y = 9) +
  scale_color_manual(values = sc_col) +
  scale_x_continuous(limits = c(-35, 35)) +
  scale_y_continuous(limits = c(-25, 45)) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none")

dev.off()

########## by treatment ##########
png(filename = "Fig2_TSNE_tg.png", width = 1250, height = 1250, res = 300)

ggplot(pt, aes(x, y, color = treatment)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  guides(color = guide_legend(title = "Treatment",
                              title.position = "left",
                              nrow = 1,
                              override.aes = list(size = 2.5,
                                                  shape = 15))) +
  scale_x_continuous(limits = c(-35, 35)) +
  scale_y_continuous(limits = c(-25, 45)) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = c(0.35, 0.9),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))

dev.off()

########## by main.class ##########
png(filename = "SFig3_TSNE_cd8.png", width = 400, height = 400, res = 300)
ggplot(pt, aes(x, y)) +
  geom_point(size = 0.1, aes(color = grepl("cd8", subclass))) +
  geom_text_repel(data = data.frame(x = -25, y = 17),
                  label = "CD8",
                  size = 2.25,
                  min.segment.length = 0,
                  nudge_x = -5,
                  nudge_y = 5) +
  scale_color_manual(values = c("grey80", "#FD5A1E")) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))
dev.off()

png(filename = "SFig3_TSNE_cd4.png", width = 400, height = 400, res = 300)
ggplot(pt, aes(x, y)) +
  geom_point(size = 0.1, aes(color = grepl("cd4", subclass))) +
  geom_text_repel(data = data.frame(x = 20, y = -5),
                  label = "CD4",
                  size = 2.25,
                  min.segment.length = 0,
                  nudge_x = 5,
                  nudge_y = -5) +
  scale_color_manual(values = c("grey80", "#AE8F6F")) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))
dev.off()

## subclustering of tnk
pt <- cbind(x = reducedDim(tnk, "UMAP")[, 1],
            y = reducedDim(tnk, "UMAP")[, 2],
            data.frame(t(logcounts(tnk)[c("CD3D", "CD8A", "CD4", "KLRD1"), ])),
            data.frame(colData(tnk)))

png(filename = "SFig2_UMAP_sc.png", width = 1000, height = 2000, res = 300)
grid.arrange(grobs = list(
  ggplot(pt, aes(x, y, color = subclass)) +
    scale_color_manual(values = sc_col) +
    geom_point(size = 0.1) +
    guides(color = guide_legend(override.aes = list(size = 1),
                                ncol = 2,
                                title = "hclust")) +
    theme_void() +
    theme(aspect.ratio = 1,
          legend.position = "none",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7))
  ,
  ggplot(pt, aes(x, y, color = treatment)) +
    geom_point(size = 0.1) +
    scale_color_manual(values = c("#1e88e5", "#d81b60")) +
    guides(color = guide_legend(override.aes = list(size = 1),
                                nrow = 2)) +
    theme_void() +
    theme(aspect.ratio = 1,
          legend.position = c(0.8, 0.8),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 7))
),
nrow = 2)
dev.off()

ts_plot <- function(gene){
  ggplot(pt, aes(x, y)) +
    geom_point(data = pt[pt[, gene] == 0, ],
               color = "grey90",
               size = 0.1) +
    geom_point(data = pt[pt[, gene] > 0, ],
               aes_string(color = gene),
               size = 0.1) +
    scale_color_gradient(low = "grey90", high = "#003831") +
    ggtitle(gene) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "italic",
                                    size = 7,
                                    hjust = 0.5))
}

png(filename = "SFig2_sm.png", width = 1000, height = 1000, res = 300)
grid.arrange(grobs = lapply(
  c("CD3D", "CD8A", "CD4", "KLRD1"),
  ts_plot)
  ,
  nrow = 2)
dev.off()

########## plot of subclass by treatment ##########
png(filename = "Fig2_sc_tg.png", width = 425, height = 1100, res = 300)

ggplot(data.frame(prop.table(table(pt$subclass, pt$treatment), 2)),
       aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "black") +
  scale_fill_distiller(palette = 7, direction = 1) +
  guides(fill = guide_colorbar(title = "Fraction of\ntreatment group",
                               title.position = "top",
                               barheight = 0.5,
                               barwidth = 3.5,
                               frame.colour = "black",
                               ticks.colour = "black")) +
  scale_x_discrete(expand = c(0,0),
                   labels = c("cont.", "SBRT")) +
  scale_y_discrete(limits = unique(pt$subclass)[order(unique(pt$subclass), decreasing = T)],
                   labels = c("NKT_GZMK", "NKT_GZMH", "NKT_KLRG1",
                              "NK_FCGR3A", "NK_NCAM1",
                              "CD8_MKI67", "CD8_HAVCR2", "CD8_CCR5", "CD8_ITM2C", "CD8_TCF7",
                              "CD4_FOXP3", "CD4_ANXA1", "CD4_FOS", "CD4_MALAT1")) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 7,
                                   color = "black"),
        axis.text.y = element_text(size = 7,
                                   color = "black",
                                   vjust = 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = "bottom",
        panel.grid = element_blank())

dev.off()

########## small multiples ##########
pt <- cbind(x = reducedDim(tnk, "TSNE")[, 1],
            y = reducedDim(tnk, "TSNE")[, 2],
            data.frame(t(logcounts(tnk)[c("CD3D", "CD8A", "CD4", "KLRD1",
                                          "TCF7", "IL7R", "HAVCR2", "TOX",
                                          "FOXP3", "IL2RA", "IFNG", "MKI67"), ])),
            data.frame(colData(tnk)))

ts_plot <- function(gene){
  ggplot(pt, aes(x, y)) +
    geom_point(data = pt[pt[, gene] == 0, ],
               color = "grey90",
               size = 0.1) +
    geom_point(data = pt[pt[, gene] > 0, ],
               aes_string(color = gene),
               size = 0.1) +
    scale_color_gradient(low = "grey90", high = "#003831") +
    ggtitle(gene) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "italic",
                                    size = 7,
                                    hjust = 0.5))
}

png(filename = "Fig2_sm.png", width = 1750, height = 1100, res = 300)

grid.arrange(grobs = lapply(
  c("CD3D", "CD8A", "CD4", "KLRD1",
    "TCF7", "IL7R", "HAVCR2", "TOX",
    "FOXP3", "IL2RA", "IFNG", "MKI67"),
  ts_plot)
,
nrow = 3)

dev.off()

## heatmap
### gene sets
s_m <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A",
         "CD8B", "NKG7", "KLRB1", "KLRD1", "KLRG1") # surface markers
n_k <- c("NCAM1", "FCGR3A") # natural killer
nai <- c("IL7R", "TCF7", "CCR7", "SELL") # naive
act <- c("CD44", "CD69", "ITM2C", "ZFP36", "TNFRSF9") # activation
exh <- c("TNFAIP3", "CTLA4", "TIGIT", "LAG3", "HAVCR2",
         "TOX", "PDCD1") # exhaustion
eff <- c("GZMK", "GZMH", "GZMB", "GZMM", "PRF1",
         "GNLY", "IFNG", "TNFSF10") # effector
trg <- c("FOXP3", "IL2RA") # t helper/treg
pro <- c("MKI67", "TOP2A")
c_e <- c("MALAT1", "NEAT1", "FOSB", "FOS", "JUNB",
         "JUN", "ANXA1", "ANXA2", "CCR5", "CCL3", "CXCL13") # cluster enriched

### subset genes
sub <- tnk[c(s_m, n_k, nai, act, exh, eff, trg, pro, c_e), ]

### calculate mean logcounts for genes of interest by cluster
heat_counts <- c()
for(x in sort(unique(tnk$subclass))){
  a <- logcounts(sub)[, which(sub$subclass == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- c("CD4_MALAT1", "CD4_FOS", "CD4_ANXA1", "CD4_FOXP3",
                           "CD8_TCF7", "CD8_ITM2C", "CD8_CCR5", "CD8_HAVCR2", "CD8_MKI67",
                           "NK_NCAM1", "NK_FCGR3A",
                           "NKT_KLRG1", "NKT_GZMH", "NKT_GZMK")

### row annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             col = list(subclass = c("CD4_MALAT1" = sc_col[1],
                                                     "CD4_FOS" = sc_col[2],
                                                     "CD4_ANXA1" = sc_col[3],
                                                     "CD4_FOXP3" = sc_col[4],
                                                     "CD8_TCF7" = sc_col[5],
                                                     "CD8_ITM2C" = sc_col[6],
                                                     "CD8_CCR5" = sc_col[7],
                                                     "CD8_HAVCR2" = sc_col[8],
                                                     "CD8_MKI67" = sc_col[9],
                                                     "NK_NCAM1" = sc_col[10],
                                                     "NK_FCGR3A" = sc_col[11],
                                                     "NKT_KLRG1" = sc_col[12],
                                                     "NKT_GZMH" = sc_col[13],
                                                     "NKT_GZMK" = sc_col[14]
                                                     )),
                             which = "column",
                             height = unit(4, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

### col annotation relative abundance
set.seed(415); cont <- colnames(tnk)[which(tnk$treatment == "control")][sample(1:8801, 5000)]
set.seed(415); sbrt <- colnames(tnk)[which(tnk$treatment == "SBRT")][sample(1:15432, 5000)]

tnk_rel <- tnk[, c(cont, sbrt)]
tnk_rel <- round(prop.table(table(tnk_rel$subclass, tnk_rel$treatment), 1) * 100, 1) %>%
  data.frame() %>%
  spread(Var2, Freq)
tnk_rel <- HeatmapAnnotation("subclass\ncomposition" = anno_barplot(cbind(tnk_rel$control, tnk_rel$SBRT),
                                                                    gp = gpar(fill = c("#1e88e5",
                                                                                       "#d81b60")),
                                                                    height = unit(1, "cm")),
                             which = "column",
                             annotation_name_gp = gpar(fontsize = 7),
                             annotation_name_rot = 0)

### heatmap scale
col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                colorRampPalette(c("purple", "black", "yellow"))(5))

########## plot heatmap ##########
png(filename = "Fig2_hm.png", width = 750, height = 2500, res = 300)

set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       top_annotation = col_ann,
                       bottom_annotation = tnk_rel,
                       show_heatmap_legend = F,
                       cluster_rows = F,
                       cluster_columns = F,
                       show_row_dend = F,
                       show_row_names = T,
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title_gp = gpar(fontsize = 7),
                       show_column_dend = F,
                       show_column_names = T,
                       column_names_side = "top",
                       column_names_gp = gpar(fontsize = 7),
                       column_title_gp = gpar(fontsize = 0),
                       column_split = factor(c(rep("CD4", 4),
                                            rep("CD8", 5),
                                            rep("NK", 2),
                                            rep("NKT", 3)),
                                          levels = c("CD4", "CD8", "NK", "NKT")),
                       row_split = factor(c(rep("lineage", 12),
                                               rep("naive", 4),
                                               rep("activation", 5),
                                               rep("exhausted", 7),
                                               rep("effector", 8),
                                               rep("Treg", 2),
                                               rep("cyc", 2),
                                               rep("cluster specific", 11)),
                                             levels = c("lineage",
                                                        "naive",
                                                        "activation",
                                                        "exhausted",
                                                        "effector",
                                                        "Treg",
                                                        "cyc",
                                                        "cluster specific")),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))

dev.off()

h_lgd = Legend(col_fun = col_fun, title = "Z score", at = c(-2, -1, 0, 1, 2),
               labels = c("< -2", "-1", "0", "1", "> 2"),
               direction = "horizontal",
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               legend_width = unit(3, "cm"),
               grid_height = unit(0.25, "cm"))
b_lgd = Legend(labels = c("control", "SBRT"),
               title = "Treatment",
               legend_gp = gpar(fill = c("#1e88e5", "#d81b60")),
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               grid_height = unit(0.25, "cm"),
               grid_width = unit(0.25, "cm"),
               nrow = 1)
lgd <- packLegend(h_lgd, b_lgd, direction = "vertical")

png(filename = "Fig2_hm_lgd.png", width = 750, height = 250, res = 300)

draw(lgd)

dev.off()

### calculate mean logcounts for genes of interest by cluster separated by treatment group
sub$tc <- paste0(sub$treatment, " ", sub$subclass, " ", "hc:", sub$hc_14)

heat_counts <- c()
for(x in sort(unique(sub$tc))){
  a <- logcounts(sub)[, which(sub$tc == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- paste(
  c(rep("control", 14), rep("SBRT", 14)),
  rep(c("CD4_MALAT1 hc:2", "CD4_FOS hc:6", "CD4_ANXA1 hc:1", "CD4_FOXP3 hc:7",
        "CD8_TCF7 hc:3", "CD8_ITM2C hc:10", "CD8_CCR5 hc:11", "CD8_HAVCR2 hc:9", "CD8_MKI67 hc:4",
        "NK_NCAM1 hc:5", "NK_FCGR3A hc:8",
        "NKT_KLRG1 hc:13", "NKT_GZMH hc:14", "NKT_GZMK hc:12"), 2))

### row annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             treatment = c(rep("control", 14), rep("SBRT",14)),
                             col = list(subclass = c("control CD4_MALAT1 hc:2" = sc_col[1],
                                                     "control CD4_FOS hc:6" = sc_col[2],
                                                     "control CD4_ANXA1 hc:1" = sc_col[3],
                                                     "control CD4_FOXP3 hc:7" = sc_col[4],
                                                     "control CD8_TCF7 hc:3" = sc_col[5],
                                                     "control CD8_ITM2C hc:10" = sc_col[6],
                                                     "control CD8_CCR5 hc:11" = sc_col[7],
                                                     "control CD8_HAVCR2 hc:9" = sc_col[8],
                                                     "control CD8_MKI67 hc:4" = sc_col[9],
                                                     "control NK_NCAM1 hc:5" = sc_col[10],
                                                     "control NK_FCGR3A hc:8" = sc_col[11],
                                                     "control NKT_KLRG1 hc:13" = sc_col[12],
                                                     "control NKT_GZMH hc:14" = sc_col[13],
                                                     "control NKT_GZMK hc:12" = sc_col[14],
                                                     "SBRT CD4_MALAT1 hc:2" = sc_col[1],
                                                     "SBRT CD4_FOS hc:6" = sc_col[2],
                                                     "SBRT CD4_ANXA1 hc:1" = sc_col[3],
                                                     "SBRT CD4_FOXP3 hc:7" = sc_col[4],
                                                     "SBRT CD8_TCF7 hc:3" = sc_col[5],
                                                     "SBRT CD8_ITM2C hc:10" = sc_col[6],
                                                     "SBRT CD8_CCR5 hc:11" = sc_col[7],
                                                     "SBRT CD8_HAVCR2 hc:9" = sc_col[8],
                                                     "SBRT CD8_MKI67 hc:4" = sc_col[9],
                                                     "SBRT NK_NCAM1 hc:5" = sc_col[10],
                                                     "SBRT NK_FCGR3A hc:8" = sc_col[11],
                                                     "SBRT NKT_KLRG1 hc:13" = sc_col[12],
                                                     "SBRT NKT_GZMH hc:14" = sc_col[13],
                                                     "SBRT NKT_GZMK hc:12" = sc_col[14]),
                                        treatment = c("control" = "#1e88e5",
                                                      "SBRT" = "#d81b60")),
                             which = "column",
                             height = unit(3, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

########## plot heatmap ##########
png(filename = "SFig2_hm.png", width = 1250, height = 2500, res = 300)
set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       top_annotation = col_ann,
                       show_heatmap_legend = F,
                       cluster_rows = F,
                       cluster_columns = F,
                       show_row_dend = F,
                       show_row_names = T,
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title_gp = gpar(fontsize = 7),
                       show_column_dend = F,
                       show_column_names = T,
                       rect_gp = gpar(col = "white", lwd = 1),
                       column_names_side = "top",
                       column_names_gp = gpar(fontsize = 7),
                       column_title_gp = gpar(fontsize = 0),
                       column_split = factor(c(rep("control", 14),
                                               rep("SBRT", 14)),
                                             levels = c("control", "SBRT")),
                       row_split = factor(c(rep("lineage", 12),
                                            rep("naive", 4),
                                            rep("activation", 5),
                                            rep("exhausted", 7),
                                            rep("effector", 8),
                                            rep("Treg", 2),
                                            rep("cyc", 2),
                                            rep("cluster specific", 11)),
                                          levels = c("lineage",
                                                     "naive",
                                                     "activation",
                                                     "exhausted",
                                                     "effector",
                                                     "Treg",
                                                     "cyc",
                                                     "cluster specific")),
                       row_gap = unit(2, "mm"),
                       column_gap = unit(5, "mm"),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))
dev.off()

## clean up
rm(list = ls()[!ls() %in% c("tnk", "sc_col")])

##########~~~CD8 t cell analysis~~~##########

## dge by limma
cd8 <- tnk[, grepl("cd8", tnk$subclass)]
cd8 <- cd8[rowSums(counts(cd8)) > 0, ] # only genes detectable within this subset
counts <- log2(1 + calculateCPM(cd8))
design <- model.matrix(~I(treatment == "SBRT"), colData(cd8))
l <- lmFit(counts, design)
e <- eBayes(l)
deg <- topTable(e, n = Inf)
deg$sig <- ifelse(abs(deg$logFC) < 1 |
                    deg$adj.P.Val > 0.05,
                  "ns",
                  ifelse(deg$logFC > 1,
                         "SBRT",
                         "cont."))

# clean up
rm(counts, design, l, e)

########## CD8 volcano plot ##########
deg_ <- deg
deg_$adj.P.Val[deg_$adj.P.Val < 10^-300] <- 10^-299.5 # rescale for cosmetics

png(filename = "Fig3_vp.png", width = 700, height = 850, res = 300)

ggplot(deg_, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.5, aes(color = sig)) +
  scale_color_manual(breaks = c("cont.", "ns", "SBRT"),
                     values = c("#1e88e5", "grey", "#d81b60")) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  shape = 15),
                              nrow = 2)) +
  labs(color = "DGE") +
  xlab(expression("control    " %<-% "    log2FC    " %->% "    SBRT")) +
  geom_text_repel(data = filter(deg_, ID %in% c("IL7R")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = 2,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg_, ID %in% c("KLF2")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = -1,
                  nudge_y = 50) +
  geom_text_repel(data = filter(deg_, ID %in% c("CCL3")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = -3,
                  nudge_y = -50) +
  geom_text_repel(data = filter(deg_, ID %in% c("JUN")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = -2.5,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg_, ID %in% c("IFNG")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = 1,
                  nudge_y = -50) +
  geom_text_repel(data = filter(deg_, ID %in% c("PDCD1", "STAT1")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = 3,
                  nudge_y = 25) +
  geom_text_repel(data = filter(deg_, ID %in% c("TOX")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = 2,
                  nudge_y = 14) +
  geom_text_repel(data = filter(deg_, ID %in% c("MKI67")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = -1,
                  nudge_y = 100) +
  geom_text_repel(data = filter(deg_, ID %in% c("CD69", "TNFSF10")),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  nudge_x = 3,
                  nudge_y = -10) +
  scale_x_continuous(limits = c(-5.5, 5.5),
                     breaks = seq(-6, 6, 2)) +
  scale_y_continuous(breaks = seq(0, 300, 100),
                     labels = c("0", "100", "200", ">300")) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0, "cm"))

dev.off()

## gsea
### reorder by logFC
g_ord <- deg[order(deg$logFC, decreasing = T), ]

### generate ENTREZID key
key <- AnnotationDbi::select(hs,
                             keys = g_ord$ID,
                             columns = c("ENTREZID"),
                             keytype = "SYMBOL")
key <- key[!duplicated(key$SYMBOL), ]

### assign ENTREZID to gene
g_ord <- left_join(g_ord, key, by = c("ID" = "SYMBOL"))

### filter duplicated and missing ENTREZID
g_ord <- g_ord[!is.na(g_ord$ENTREZID), ]
g_ord <- g_ord[!duplicated(g_ord$ENTREZID), ]

### generate ENTREZID list
ent_ord <- g_ord[, "logFC"]
names(ent_ord) <- as.character(g_ord[, "ENTREZID"])

### reactome pathway enrichment
set.seed(415); reactomeGSEA <- gsePathway(ent_ord,
                                          maxGSSize = 500,
                                          pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH",
                                          verbose = FALSE) %>%
  data.frame()

sub_GSEA <- reactomeGSEA %>%
  filter(grepl("Mitotic|Interferon|interferon|TCR", Description))
sub_GSEA$Description <- factor(sub_GSEA$Description,
                               levels = sub_GSEA$Description[order(sub_GSEA$NES, decreasing = F)])

########## CD8 GSEA ##########
png(filename = "Fig3_gsea.png", width = 600, height = 800, res = 300)

ggplot(sub_GSEA, aes(x = Description, y = -log10(p.adjust))) +
  geom_col(aes(alpha = NES),
           fill = "#d81b60",
           color = "darkgrey") +
  guides(alpha = guide_legend(nrow = 2)) +
  geom_text(aes(label = Description),
            y = 0.025,
            hjust = "inward",
            size = 2.25) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "twodash",
             alpha = 0.50) +
  xlab("Reactome Pathway Enrichment") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.text.x = element_text(size = 7, color = "black"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.25, 'cm'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")

dev.off()

## pseudotime
### recalculate UMAP for CD8 only
set.seed(415); cd8 <- runUMAP(cd8)

### calculate UMAP, use existing subcluster id
cd8 <- slingshot(cd8, clusterLabels = "subclass", reducedDim = "UMAP")

cd8$Pseudotime <- ifelse(is.na(cd8$slingPseudotime_1),
                         cd8$slingPseudotime_2,
                         cd8$slingPseudotime_1)

pt <- cbind(x = -reducedDims(cd8)$UMAP[, 1],
            y = reducedDims(cd8)$UMAP[, 2],
            data.frame(colData(cd8)[, c("treatment", "subclass", "sample",
                                        "Pseudotime", "slingPseudotime_1", "slingPseudotime_2")]),
            data.frame(t(logcounts(cd8)[c("TCF7", "CCR7", "IL7R",
                                          "CD69", "ZFP36", "ITM2C",
                                          "PDCD1", "HAVCR2", "LAG3",
                                          "IFNG", "TNFSF10", "GZMK",
                                          "TOX", "MKI67", "TOP2A",
                                          "CCL3", "CCR5", "ANXA2",
                                          "JUN", "STAT1", "PRF1"), ])))

### calculate lineage
lin <- getLineages(reducedDims(cd8)$UMAP, cd8$subclass)
cur <- getCurves(lin)

cur1 <- data.frame(x = -cur@metadata$curves$Lineage1$s[, 1],
                   y = cur@metadata$curves$Lineage1$s[, 2])
cur2 <- data.frame(x = -cur@metadata$curves$Lineage2$s[, 1],
                   y = cur@metadata$curves$Lineage2$s[, 2])

# clean up
rm(lin, cur)

########## CD8 pseudotime ##########
png(filename = "Fig3_ps_sc.png", width = 1000, height = 450, res = 300)
ggplot(pt, aes(x, y)) +
  scale_color_manual(values = sc_col[5:9]) +
  geom_point(aes(color = subclass),
             size = 0.25) +
  geom_point(data = cur1,
             size = 0.1,
             shape = 46) +
  geom_point(data = cur2,
             size = 0.1,
             shape = 46) +
  geom_text_repel(data = data.frame(x = -3.5, y = 2.75),
                  label = "CD8_TCF7",
                  min.segment.length = 0,
                  nudge_x = -0.25,
                  nudge_y = -2,
                  size = 2.25) +
  geom_text_repel(data = data.frame(x = -3, y = -0.75),
                  label = "CD8_ITM2C",
                  min.segment.length = 0,
                  nudge_x = -0.5,
                  nudge_y = -1,
                  size = 2.25) +
  geom_text_repel(data = data.frame(x = -1, y = -3.5),
                  label = "CD8_HAVCR2",
                  min.segment.length = 0,
                  nudge_x = -2,
                  nudge_y = 0,
                  size = 2.25) +
  geom_text_repel(data = data.frame(x = 1.25, y = 3.75),
                  label = "CD8_CCR5",
                  min.segment.length = 0,
                  nudge_x = 1.5,
                  nudge_y = 0,
                  size = 2.25) +
  geom_text_repel(data = data.frame(x = 5.5, y = 0.75),
                  label = "CD8_MKI67",
                  min.segment.length = 0,
                  nudge_x = 1,
                  nudge_y = 1.5,
                  size = 2.25) +
  theme_void() +
  theme(aspect.ratio = 0.33,
        title = element_text(size = 7),
        plot.title = element_text(size = 7,
                                  hjust = 0.5),
        legend.text = element_text(size = 7),
        legend.position = "none")
dev.off()

png(filename = "Fig3_ps_tg.png", width = 1450, height = 900, res = 300)

grid.arrange(grobs = list(
  ggplot(pt, aes(x, y)) +
    geom_point(aes(color = slingPseudotime_2),
               size = 0.1) +
    scale_color_distiller(palette = "BuGn", direction = 1, limits = c(0, 15), na.value = "grey75") +
    guides(color = guide_colorbar(barwidth = 5, barheight = 0.25,
                                  frame.colour = "black", ticks.colour = "black",
                                  title = "CCR5 lineage\npseudotime",
                                  title.position = "left")) +
    geom_point(data = cur1,
               size = 0.1,
               shape = 46) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          title = element_text(size = 7),
          plot.title = element_text(size = 7,
                                    hjust = 0.5),
          legend.text = element_text(size = 7),
          legend.position = "bottom")
  ,
  ggplot(pt, aes(x, y)) +
    geom_point(aes(color = slingPseudotime_1),
               size = 0.1) +
    scale_color_distiller(palette = "BuPu", direction = 1, limits = c(0, 15), na.value = "grey75") +
    guides(color = guide_colorbar(barwidth = 5, barheight = 0.25,
                                  frame.colour = "black", ticks.colour = "black",
                                  title = "MKI67 lineage\npseudotime",
                                  title.position = "left")) +
    geom_point(data = cur1,
               size = 0.1,
               shape = 46) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          title = element_text(size = 7),
          plot.title = element_text(size = 7,
                                    hjust = 0.5),
          legend.text = element_text(size = 7),
          legend.position = "bottom")
  ,
  ggplot(pt, aes(x, y)) +
    ggtitle("control\n") +
    geom_point(data = filter(pt, treatment == "SBRT"),
               color = "grey75",
               size = 0.1) +
    geom_point(data = filter(pt, treatment == "control"),
               color = "#1e88e5",
               size = 0.1) +
    geom_point(data = cur1,
               size = 0.1,
               shape = 46) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          title = element_text(size = 7),
          plot.title = element_text(size = 7,
                                    hjust = 0.5))
  ,
  ggplot(pt, aes(x, y)) +
    ggtitle("SBRT\n") +
    geom_point(data = filter(pt, treatment == "control"),
               color = "grey75",
               size = 0.1) +
    geom_point(data = filter(pt, treatment == "SBRT"),
               color = "#d81b60",
               size = 0.1) +
    geom_point(data = cur1,
               size = 0.1,
               shape = 46) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          plot.title = element_text(size = 7,
                                    hjust = 0.5))
),
layout_matrix = rbind(c(1, 2), c(3, 4)))

dev.off()

########## CD8 pseudotime density by treatment ##########
png(filename = "Fig3_ps_dn.png", width = 900, height = 600, res = 300)

grid.arrange(grobs = list(
  ggplot(filter(pt, treatment == "control"), aes(x = slingPseudotime_2)) +
    geom_density(fill = "#1e88e5") +
    ggtitle("control") +
    ylab("CCR5 lineage\ndensity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(plot.title = element_text(size = 7,
                                    color = "black",
                                    hjust = 0.5),
          axis.title.x = element_text(size = 7,
                                      color = "white"),
          axis.title.y = element_text(size = 7,
                                      color = "black"),
          axis.text.x = element_text(size = 7,
                                     color = "white"),
          axis.text.y = element_text(size = 7,
                                     color = "black"),
          panel.grid.major.x = element_line(color = "grey90"))
  ,
  ggplot(filter(pt, treatment == "SBRT"), aes(x = slingPseudotime_2)) +
    geom_density(fill = "#d81b60") +
    ggtitle("SBRT") +
    ylab("CCR5 lineage\ndensity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1)) +
    theme_classic() +
    theme(plot.title = element_text(size = 7,
                                    color = "black",
                                    hjust = 0.5),
          axis.title.x = element_text(size = 7,
                                      color = "white"),
          axis.title.y = element_text(size = 7,
                                      color = "white"),
          axis.text.x = element_text(size = 7,
                                     color = "white"),
          axis.text.y = element_text(size = 7,
                                     color = "black"),
          panel.grid.major.x = element_line(color = "grey90"))
  ,
  ggplot(filter(pt, treatment == "control"), aes(x = slingPseudotime_1)) +
    geom_density(fill = "#1e88e5") +
    ggtitle("control") +
    xlab("Pseudotime") +
    ylab("MKI67 lineage\ndensity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4)) +
    theme_classic() +
    theme(plot.title = element_text(size = 7,
                                    color = "white",
                                    hjust = 0.5),
          axis.title = element_text(size = 7,
                                    color = "black"),
          axis.text = element_text(size = 7,
                                   color = "black"),
          panel.grid.major.x = element_line(color = "grey90"))
  ,
  ggplot(filter(pt, treatment == "SBRT"), aes(x = slingPseudotime_1)) +
    geom_density(fill = "#d81b60") +
    ggtitle("SBRT") +
    xlab("Pseudotime") +
    ylab("MKI67 lineage\ndensity") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 15)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1)) +
    theme_classic() +
    theme(plot.title = element_text(size = 7,
                                    color = "white",
                                    hjust = 0.5),
          axis.title.x = element_text(size = 7,
                                      color = "black"),
          axis.title.y = element_text(size = 7,
                                      color = "white"),
          axis.text.x = element_text(size = 7,
                                     color = "black"),
          axis.text.y = element_text(size = 7,
                                     color = "black"),
          panel.grid.major.x = element_line(color = "grey90"))
)
, nrow = 2)

dev.off()

########## cd8 starting vs terminal node ##########
# subclass proportion by sample
## calculate subclass fractions by sample
samp_comp <- data.frame(subclass = unique(pt$subclass))
for(i in unique(pt$sample)){
  pt_ <- filter(pt, sample == i)
  samp_comp_ <- data.frame(prop.table(table(subclass = pt_$subclass)))
  rownames(samp_comp_) <- samp_comp_$Var1
  colnames(samp_comp_)[2] <- paste0(i)
  samp_comp <- left_join(samp_comp, samp_comp_, by = "subclass")
  rm(pt_, samp_comp_)
}
rownames(samp_comp) <- samp_comp$subclass
samp_comp <- data.frame(t(as.matrix(samp_comp[, 2:ncol(samp_comp)])))
samp_comp$sample <- rownames(samp_comp)

## rearrange for plot
samp_comp <- left_join(samp_comp, distinct(pt[, c("sample", "treatment")]), by = "sample")

### Naive node: CD8_TCF7
samp_comp_ <- samp_comp[, c("sample", "cd8_1", "treatment")]
samp_comp_ <- pivot_longer(samp_comp_, cols = "cd8_1", names_to = "subclass", values_to = "fraction")
samp_comp_$subclass <- factor(samp_comp_$subclass, levels = c("cd8_1"))

ggplot(samp_comp_, aes(fill = subclass, y = fraction, x = sample)) +
  scale_fill_manual(values = c("#1e88e5", "#d81b60")) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = 16),
                             keyheight = unit(4, "mm"),
                             keywidth = unit(4, "mm")),
         color = guide_legend(override.aes = list(size = 1, shape = 16),
                              keyheight = unit(4, "mm"),
                              keywidth = unit(4, "mm"))) +
  geom_bar(position = "stack", stat = "identity", aes(fill = treatment)) +
  ylab("Fraction of CD8 T cells in starting node\n(CD8_TCF7)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.75)) +
  theme_classic() +
  theme(aspect.ratio = 2.5,
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        legend.position = "none") -> cd8_start

### Terminal nodes: CD8_MKI67, CD8_CCR5
samp_comp_ <- samp_comp[, c("sample", "cd8_3", "cd8_5", "treatment")]
samp_comp_ <- pivot_longer(samp_comp_, cols = c("cd8_3", "cd8_5"), names_to = "subclass", values_to = "fraction")
samp_comp_$subclass <- factor(samp_comp_$subclass, levels = c("cd8_3", "cd8_5"))

ggplot(samp_comp_, aes(fill = subclass, y = fraction, x = sample)) +
  scale_fill_manual(values = c("#1e88e5", "#d81b60")) +
  guides(fill = guide_legend(override.aes = list(size = 1, shape = 16),
                             keyheight = unit(4, "mm"),
                             keywidth = unit(4, "mm")),
         color = guide_legend(override.aes = list(size = 1, shape = 16),
                              keyheight = unit(4, "mm"),
                              keywidth = unit(4, "mm"))) +
  geom_bar(position = "stack", stat = "identity", aes(fill = treatment)) +
  ylab("Fraction of CD8 T cells in terminal nodes\n(CD8_CCR5, CD8_MKI67)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.75)) +
  theme_classic() +
  theme(aspect.ratio = 2.5,
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        legend.position = "none") -> cd8_term

png("cd8_term_clustfrac.png", height = 500, width = 500, res = 150)
grid.arrange(grobs = list(cd8_start, cd8_term), nrow = 1)
dev.off()

########## cd8 imm/act chisq ##########
set.seed(415); cont <- colnames(cd8)[which(cd8$treatment == "control")][sample(1:1528, 1250)]
set.seed(415); sbrt <- colnames(cd8)[which(cd8$treatment == "SBRT")][sample(1:7335, 1250)]
cd8_1250 <- cd8[, c(cont, sbrt)]

chisq.test(table(cd8_1250$treatment, cd8_1250$subclass == "cd8_1"))
cs <- data.frame(table(treatment = cd8_1250$treatment, naive = cd8_1250$subclass == "cd8_1"))
cs$x0 <- c(2, 2, 1, 1)
cs$y0 <- c(2, 1, 2, 1)
cs$frac <- cs$Freq/max(cs$Freq)
cs$Freq <- cs$Freq/sum(cs$Freq)

png(filename = "Fig5_ps_cs.png", width = 480, height = 600, res = 300)

ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = frac/2, fill = Freq), data = cs) +
  scale_fill_distiller(palette = 7, direction = 1) +
  guides(fill = guide_colorbar(title = "Fraction of cells",
                               title.position = "top",
                               barheight = 0.5,
                               barwidth = 3.5,
                               frame.colour = "black",
                               ticks.colour = "black",
                               label.theme = element_text(angle = 90,
                                                          vjust = 0.5,
                                                          size = 7))) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("naive", "activated")) +
  scale_y_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("SBRT", "control")) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_text(size = 7,
                                 color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        legend.title = element_text(size = 7))

dev.off()

########## CD8 pseudotime smoothed expression by lineage ##########

pt_mat <- gather(pt,
                 "IL7R", "TCF7", "CCR7",
                 "PDCD1", "TOX", "HAVCR2",
                 "CCR5", "MKI67", "TOP2A",
                 key = "gene", value = "logcounts")
pt_mat$gene <- factor(pt_mat$gene,
                      levels = c("IL7R", "TCF7", "CCR7",
                                 "PDCD1", "TOX", "HAVCR2",
                                 "CCR5", "MKI67", "TOP2A"))

png(filename = "Fig3_ps_sm_ln_mat.png", width = 1450, height = 700, res = 300)

ggplot(pt_mat,
       aes(y = logcounts)) +
  geom_smooth(aes(x = slingPseudotime_2), se = F, color = "seagreen3") +
  geom_smooth(aes(x = slingPseudotime_1), se = F, color = "darkorchid3") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Pseudotime") +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7, face = "italic",
                                  margin = margin(t = 1, b = 1)),
        strip.background = element_rect(fill = "white")) +
  theme(legend.position = "right",
        aspect.ratio = 0.33) +
  facet_wrap(~gene,
             scale = "free_y")

dev.off()

# clean up
rm(pt_mat)

########## CD8 pseudotime smoothed expression by treatment ##########

pt_tg <- gather(pt,
                "CD69", "JUN", "STAT1",
                "TNFSF10", "IFNG", "PRF1",
                 key = "gene", value = "logcounts")
pt_tg$gene <- factor(pt_tg$gene,
                     levels = c("CD69", "JUN", "STAT1",
                                "TNFSF10", "IFNG", "PRF1"))

png(filename = "Fig3_ps_sm_tg.png", width = 1450, height = 500, res = 300)

ggplot(pt_tg,
       aes(y = logcounts)) +
  geom_smooth(aes(x = Pseudotime, color = treatment), se = F) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Pseudotime") +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(size = 7, face = "italic",
                                  margin = margin(t = 1, b = 1)),
        strip.background = element_rect(fill = "white")) +
  theme(legend.position = "none",
        aspect.ratio = 0.33) +
  facet_wrap(~gene,
             scale = "free_y")

dev.off()

l_lgd <- Legend(labels = c("control", "SBRT"),
               title = "Treatment",
               legend_gp = gpar(fill = c("#1e88e5", "#d81b60")),
               title_gp = gpar(fontsize = 7),
               title_position = "leftcenter",
               labels_gp = gpar(fontsize = 7),
               grid_height = unit(0.25, "cm"),
               grid_width = unit(0.25, "cm"),
               nrow = 1)
t_lgd <- Legend(labels = c("CCR5", "MKI67"),
                title = "Lineage",
                legend_gp = gpar(fill = c("seagreen3", "darkorchid3")),
                title_gp = gpar(fontsize = 7),
                title_position = "leftcenter",
                labels_gp = gpar(fontsize = 7),
                grid_height = unit(0.25, "cm"),
                grid_width = unit(0.25, "cm"),
                nrow = 1)
lgd <- packLegend(l_lgd, t_lgd, direction = "vertical", row_gap = unit(1, "cm"))

png(filename = "Fig2_hm_lgd.png", width = 1450, height = 250, res = 300)

draw(lgd)

dev.off()

########## CD8 gene expression by UMAP ##########
ps_plot <- function(gene){
  ggplot(pt, aes(x, y)) +
    geom_point(data = pt[pt[, gene] == 0, ],
               color = "grey90",
               size = 0.1) +
    geom_point(data = pt[pt[, gene] > 0, ],
               aes_string(color = gene,
                          alpha = gene),
               size = 0.1) +
    scale_color_gradient2(low = "grey90", mid = "#003831", high = "#efb21e") +
    guides(color = guide_colorbar(barwidth = 3,
                                  barheight = 0.25,
                                  title.position = "left"),
           alpha = F) +
    geom_point(data = cur1,
               size = 0.1,
               shape = 46) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          legend.title = element_text(size = 7, hjust = 0, face = "italic"),
          legend.text = element_text(size = 7),
          legend.position = "bottom")
}

png(filename = "Fig3_ps_tg.png", width = 1450, height = 1800, res = 300)

grid.arrange(grobs = lapply(
  c("TCF7", "CCR7", "IL7R",
    "CD69", "ZFP36", "ITM2C",
    "PDCD1", "HAVCR2", "LAG3",
    "IFNG", "TNFSF10", "GZMK",
    "TOX", "MKI67", "TOP2A",
    "CCL3", "CCR5", "ANXA2"),
  ps_plot
)
,
nrow = 6)

dev.off()

# clean up
rm(list = ls()[!ls() %in% c("tnk", "hs")])

##########~~~CD4 t cell analysis~~~##########

## dge by limma
cd4 <- tnk[, grepl("cd4", tnk$subclass)]
cd4 <- cd4[rowSums(counts(cd4)) > 0, ] # only genes detectable within this subset
counts <- log2(1 + calculateCPM(cd4))
design <- model.matrix(~I(treatment == "SBRT"), colData(cd4))
l <- lmFit(counts, design)
e <- eBayes(l)
deg <- topTable(e, n = Inf)
deg$sig <- ifelse(abs(deg$logFC) < 1 |
                    deg$adj.P.Val > 0.05,
                  "ns",
                  ifelse(deg$logFC > 1,
                         "SBRT",
                         "cont."))

# clean up
rm(counts, design, l, e)

########## CD4 volcano plot ##########
png(filename = "SFig3_vp.png", width = 1000, height = 1450, res = 300)

ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.5, aes(color = sig)) +
  scale_color_manual(breaks = c("cont.", "ns", "SBRT"),
                     values = c("#1e88e5", "grey", "#d81b60")) +
  scale_x_continuous(limits = c(-4, 4)) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  shape = 15))) +
  geom_text_repel(data = filter(top_n(deg, -10, adj.P.Val)),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  min.segment.length = 0,
                  xlim = c(2.5, NA),
                  direction = "y") +
  geom_text_repel(data = filter(top_n(deg, -10, logFC)),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  min.segment.length = 0,
                  xlim = c(NA, -2.5),
                  direction = "y") +
  geom_text_repel(data = filter(deg, grepl("FOXP3|STAT1|STAT3", ID)),
                  aes(label = ID), fontface = "italic",
                  size = 2.25,
                  min.segment.length = 0,
                  xlim = c(2.5, NA),
                  direction = "y") +
  theme_classic() +
  xlab(expression("control    " %<-% "    log2FC    " %->% "    SBRT")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0, "cm"))

dev.off()

## gsea
### reorder by logFC
g_ord <- deg[order(deg$logFC, decreasing = T), ]

### generate ENTREZID key
key <- AnnotationDbi::select(hs,
                             keys = g_ord$ID,
                             columns = c("ENTREZID"),
                             keytype = "SYMBOL")
key <- key[!duplicated(key$SYMBOL), ]

### assign ENTREZID to gene
g_ord <- left_join(g_ord, key, by = c("ID" = "SYMBOL"))

### filter duplicated and missing ENTREZID
g_ord <- g_ord[!is.na(g_ord$ENTREZID), ]
g_ord <- g_ord[!duplicated(g_ord$ENTREZID), ]

### generate ENTREZID list
ent_ord <- g_ord[, "logFC"]
names(ent_ord) <- as.character(g_ord[, "ENTREZID"])

### kegg pathway enrichment
set.seed(415); keggGSEA <- gseKEGG(ent_ord,
                                   maxGSSize = 500,
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH",
                                   verbose = FALSE) %>%
  data.frame()

sub_GSEA <- keggGSEA %>%
  filter(grepl("Th1|PD", Description))
sub_GSEA$Description <- factor(sub_GSEA$Description,
                               levels = sub_GSEA$Description[order(sub_GSEA$NES, decreasing = F)])

########## CD4 GSEA ##########
png(filename = "CD4_gsea.png", width = 850, height = 400, res = 300)

ggplot(sub_GSEA, aes(x = Description, y = -log10(p.adjust))) +
  geom_col(aes(alpha = NES),
           fill = "#d81b60",
           color = "darkgrey") +
  guides(alpha = guide_legend(nrow = 2)) +
  geom_text(aes(label = Description),
            y = 0.025,
            hjust = "inward",
            size = 2.25) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "twodash",
             alpha = 0.50) +
  ggtitle("KEGG Pathway Enrichment") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.25, 'cm'),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")

dev.off()

########## cd4 ANXA1 chisq ##########
chisq.test(table(cd4$treatment, cd4$subclass == "cd4_3"))
cs <- data.frame(table(treatment = cd4$treatment, CD4_ANXA1 = cd4$subclass == "cd4_3"))
cs$x0 <- c(2, 2, 1, 1)
cs$y0 <- c(2, 1, 2, 1)
cs$frac <- cs$Freq/max(cs$Freq)
cs$Freq <- cs$Freq/sum(cs$Freq)

png(filename = "SFig2_ps_cs.png", width = 480, height = 600, res = 300)

ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = frac/2, fill = Freq), data = cs) +
  scale_fill_distiller(palette = 7, direction = 1) +
  guides(fill = guide_colorbar(title = "Fraction of cells",
                               title.position = "top",
                               barheight = 0.5,
                               barwidth = 3.5,
                               frame.colour = "black",
                               ticks.colour = "black",
                               label.theme = element_text(angle = 90,
                                                          vjust = 0.5,
                                                          size = 7))) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("naive", "activated")) +
  scale_y_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("SBRT", "control")) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_text(size = 7,
                                 color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.5),
        legend.title = element_text(size = 7))

dev.off()
