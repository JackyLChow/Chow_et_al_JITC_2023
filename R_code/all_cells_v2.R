# load packages
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(viridis)
library(readr)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(tidyr)
library(kableExtra)

# load data
sce <- readRDS("~/_Projects/Chow Nat Comm/data/sce_classed_FZ02.RDS")
a_c <- cbind(data.frame(reducedDim(sce, "TSNE")),
             colData(sce)[, c("treatment", "class", "sample")])
a_c$subclass <- NA
a_c$subclass_fine <- NA

# assign subclasses from tnk
tnk <- readRDS("~/_Projects/Chow Nat Comm/data/tnk_FZ02.RDS")
a_c[colnames(tnk), ]$subclass <- tnk$subclass
a_c[colnames(tnk), ]$subclass_fine <- tnk$subclass_fine

# cell x sample table
#write.csv(table(tnk$subclass_fine, tnk$sample), "TableS2_cellxsample_tnk.csv")

rm(tnk) # clean up

# assign subclasses from mye
mye <- readRDS("~/_Projects/Chow Nat Comm/data/mye_FZ02.RDS")
a_c[colnames(mye), ]$subclass <- mye$subclass
a_c[colnames(mye), ]$subclass_fine <- mye$subclass_fine

# cell x sample table
#write.csv(table(mye$subclass_fine, mye$sample), "TableS3_cellxsample_mye.csv")

rm(mye) # clean up

# inherit b_c, ves, tum, and mst from class
a_c[a_c$class %in% c("b_c", "ves", "tum", "mst"), ]$subclass <- a_c[a_c$class %in% c("b_c", "ves", "tum", "mst"), ]$class
a_c[a_c$class %in% c("b_c", "ves", "tum", "mst"), ]$subclass_fine <- a_c[a_c$class %in% c("b_c", "ves", "tum", "mst"), ]$class

a_c$subclass_main <- substr(a_c$subclass, 1, 3)

# add annotation to sce
sce$subclass <- a_c$subclass
sce$subclass_main <- a_c$subclass_main
sce$subclass_fine <- a_c$subclass_fine

# save
#saveRDS(sce, "~/_Projects/Chow Nat Comm/data/sce_final.RDS")

# break up for upload
for(brk in 1:ceiling(ncol(sce)/8000)){
  brk_min_ <- 1 + (8000 * (brk - 1))
  brk_max_ <- 8000 + (8000 * (brk - 1))
  if(brk_max_ < ncol(sce)){
    sce_ <- sce[, brk_min_:brk_max_]
    saveRDS(sce_, paste0("~/_Projects/Chow Nat Comm/data/sce_final/sce_final_chunk_", brk, ".RDS"))
  } else {
    sce_ <- sce[, brk_min_:ncol(sce)]
    saveRDS(sce_, paste0("~/_Projects/Chow Nat Comm/data/sce_final/sce_final_chunk_", brk, ".RDS"))
  }
  rm(sce_, brk) #clean up
}

# final, fully annotated sce, contains TSNE and UMAP calculated for all cells
sce <- SingleCellExperiment()
for(file_ in list.files("~/_Projects/Chow Nat Comm/data/sce_final", full.names = T)){
  if(nrow(sce) == 0){
    sce <- readRDS(file_)
  } else {
    sce <- cbind(sce, readRDS(file_))
  }
}

# cell x sample table
write.csv(table(a_c$subclass_main, a_c$sample), "TableS1_cellxsample.csv")

# genes of interest
a_c <- cbind(a_c, 
             data.frame(t(logcounts(sce)[c("PTPRC",
                                           "CD3D", "CD3E", "CD3G",
                                           "CD8A", "CD8B",
                                           "CD4",
                                           "NCAM1", "FCGR3A", "NKG7",
                                           "S100A8", "S100A9",
                                           "CD68", "HLA-DRA",
                                           "CLEC9A", "CLEC10A", "FCER1A",
                                           "CD19", "MS4A1",
                                           "CA9", "NAT8"), ])))

a_c <- data.frame(colData(sce),
                  data.frame(t(logcounts(sce)[c("PTPRC",
                                                "CD3D", "CD3E", "CD3G",
                                                "CD8A", "CD8B",
                                                "CD4",
                                                "NCAM1", "FCGR3A", "NKG7",
                                                "S100A8", "S100A9",
                                                "CD68", "HLA-DRA",
                                                "CLEC9A", "CLEC10A", "FCER1A",
                                                "CD19", "MS4A1",
                                                "CA9", "NAT8"), ])),
                  
                  reducedDim(sce, "TSNE"))

# plot TSNE and small multiples
png("Fig1_tsne.png", width = 1290, height = 1290, res = 300)
ggplot(a_c) +
  geom_point(aes(x = X1, y = X2, color = substr(subclass, 1, 3)), size = 0.2) +
  geom_label_repel(data = data.frame(x = -23, y = 18), aes(x, y),
                  color = "#cc5803",
                  segment.color = "black",
                  label = "CD8",
                  min.segment.length = 0,
                  nudge_x = -5,
                  nudge_y = 5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = -5, y = 30), aes(x, y),
                  color = "#ffb627",
                  segment.color = "black",
                  label = "NKT",
                  min.segment.length = 0,
                  nudge_x = -5,
                  nudge_y = 5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 17, y = 24), aes(x, y),
                  color = "#e2711d",
                  segment.color = "black",
                  label = "NK",
                  min.segment.length = 0,
                  nudge_x = 5,
                  nudge_y = 5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 20, y = 10), aes(x, y),
                  color = "#ff9505",
                  segment.color = "black",
                  label = "CD4",
                  min.segment.length = 0,
                  nudge_x = 5,
                  nudge_y = 5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 5, y = -16), aes(x, y),
                  color = "#3a0ca3",
                  segment.color = "black",
                  label = "DC",
                  min.segment.length = 0,
                  nudge_x = -12,
                  nudge_y = -2,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 16, y = -11), aes(x, y),
                  color = "#1b9aaa",
                  segment.color = "black",
                  label = "Mast",
                  min.segment.length = 0,
                  nudge_x = 17,
                  nudge_y = 10,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 28, y = -17), aes(x, y),
                  color = "#4cc9f0",
                  segment.color = "black",
                  label = "Monocyte",
                  min.segment.length = 0,
                  nudge_x = 5,
                  nudge_y = -5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = 5, y = -33), aes(x, y),
                  color = "#4361ee",
                  segment.color = "black",
                  label = "Macrophage",
                  min.segment.length = 0,
                  nudge_x = 5,
                  nudge_y = -5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = -14, y = -23), aes(x, y),
                  color = "#d00000",
                  segment.color = "black",
                  label = "Vessel",
                  min.segment.length = 0,
                  nudge_x = -5,
                  nudge_y = -5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = -14, y = -33), aes(x, y),
                  color = "#06d6a0",
                  segment.color = "black",
                  label = "B cell",
                  min.segment.length = 0,
                  nudge_x = -5,
                  nudge_y = -5,
                  size = 2.25) +
  geom_label_repel(data = data.frame(x = -23, y = -24), aes(x, y),
                  color = "#424c55",
                  segment.color = "black",
                  label = "Tumor",
                  min.segment.length = 0,
                  nudge_x = -7,
                  nudge_y = -5,
                  size = 2.25) +
  scale_color_manual(breaks = c("cd4", "cd8", "n_k", "nkt",
                                "mon", "mac", "d_c",
                                "mst",
                                "b_c",
                                "tum", "ves"),
                     values = c("#ff9505", "#cc5803", "#e2711d", "#ffb627",
                                "#4cc9f0", "#4361ee", "#3a0ca3",
                                "#1b9aaa",
                                "#06d6a0",
                                "#424c55", "#d00000")) +
  guides(color = guide_legend(title = "class",
                              override.aes = list(size = 3,
                                                  shape = 15),
                              nrow = 2,
                              title.position = "top")) +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
dev.off()

png("Fig1_tg.png", width = 900, height = 600, res = 300)
ggplot(a_c, aes(X1, X2, color = treatment)) +
  geom_point(size = 0.1) +
  scale_color_manual(breaks = c("control", "SBRT"),
                     values = c("#1e88e5", "#d81b60")) +
  guides(color = guide_legend(title = "",
                              override.aes = list(size = 3,
                                                  shape = 15),
                              nrow = 2,
                              title.position = "top")) +
  theme_void() +
  theme(legend.position = "right",
        aspect.ratio = 1,
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7))
dev.off()

png("Fig1_sm.png", width = 800, height = 800, res = 300)
grid.arrange(grobs = list(
ggplot(a_c, aes(X1, X2, color = PTPRC)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "grey90", high = "#003831") +
  ggtitle("PTPRC") +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(face = "italic",
                                  size = 7,
                                  hjust = 0.5))
,
ggplot(a_c, aes(X1, X2, color = CA9)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "grey90", high = "#003831") +
  ggtitle("CA9") +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(face = "italic",
                                  size = 7,
                                  hjust = 0.5))
,
ggplot(a_c, aes(X1, X2, color = CD3G)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "grey90", high = "#003831") +
  ggtitle("CD3G") +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(face = "italic",
                                  size = 7,
                                  hjust = 0.5))
,
ggplot(a_c, aes(X1, X2, color = CD68)) +
  geom_point(size = 0.1) +
  scale_color_gradient(low = "grey90", high = "#003831") +
  ggtitle("CD68") +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        plot.title = element_text(face = "italic",
                                  size = 7,
                                  hjust = 0.5))
),
nrow = 2
)
dev.off()

png("Fig1_el.png", width = 350, height = 120, res = 300)
exp_lgd <- Legend(col_fun = circlize::colorRamp2(seq(0, 5, 0.5),
                                                 colorRampPalette(c("grey90", "#003831"))(11)),
                  title = "Gene expression",
                  at = c(0, 5),
                  labels = c("low", "high"),
                  direction = "horizontal",
                  title_gp = gpar(fontsize = 7),
                  labels_gp = gpar(fontsize = 7),
                  legend_width = unit(2, "cm"),
                  grid_height = unit(0.25, "cm"))
draw(exp_lgd)
dev.off()

sce$subclass_main <- substr(a_c$subclass, 1, 3) # apply subclasses
sce$subclass.fine <- a_c$subclass 

# plot classifying heatmap
sub <- sce[c("PTPRC",
             "CD3D", "CD3E", "CD3G",
             "CD8A", "CD8B",
             "CD4", "FOXP3", "IL2RA",
             "NKG7", "KLRD1", "KLRG1", "NCAM1",
             "CD68", "HLA-DRA", 
             "S100A8", "S100A9","APOE", "APOC1",
             "CLEC9A", "CLEC10A", "FCER1A",
             "KIT", "TPSAB1",
             "CD19", "MS4A1",
             "CA9", "NAT8",
             "PECAM1", "CD34"), ]

heat_counts <- c()
for(x in sort(unique(sce$subclass_main))){
  a <- logcounts(sub)[, which(sub$subclass_main == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2
heat_counts <- heat_counts[, c("cd4","cd8", "n_k", "nkt",
                               "mon", "mac", "d_c", "mst", "b_c",
                               "tum", "ves")]
colnames(heat_counts) <- c("CD4", "CD8", "NK", "NKT",
                           "Monocyte", "Macrophage", "DC", "Mast", "B cell",
                           "Tumor", "Vessel")

col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                colorRampPalette(c("purple", "black", "yellow"))(5))

col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             col = list(subclass = c("CD4" = "#ff9505",
                                                     "CD8" = "#cc5803",
                                                     "NK" = "#e2711d",
                                                     "NKT" = "#ffb627",
                                                     "Monocyte" = "#4cc9f0",
                                                     "Macrophage" = "#4361ee",
                                                     "DC" = "#3a0ca3",
                                                     "Mast" = "#1b9aaa",
                                                     "B cell" = "#06d6a0",
                                                     "Tumor" = "#424c55",
                                                     "Vessel" = "#d00000"
                             )),
                             which = "column",
                             gp = gpar(col = "white"),
                             simple_anno_size = unit(3, "mm"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

png("Fig1_hm.png", width = 600, height = 1400, res = 300)
Heatmap(heat_counts,
        col = col_fun,
        border = F,
        top_annotation = col_ann,
        row_names_side = "right",
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = F,
        cluster_columns = F,
        width = ncol(heat_counts)*unit(3, "mm"), 
        height = nrow(heat_counts)*unit(3, "mm"),
        show_heatmap_legend = F)
dev.off()

rm(col_ann, heat_counts, col_fun, x, a)

# class.main composition
tab <- data.frame(table(subclass_main = sce$subclass_main, treatment = sce$treatment))
tab$subclass_main <- factor(tab$subclass_main,
                            levels = c("cd4", "cd8", "n_k", "nkt",
                                       "mon", "mac", "d_c",
                                       "mst",
                                       "b_c",
                                       "tum", "ves"))

png("SFig1_cc.png", width = 950, height = 1400, res = 300)
ggplot(tab, aes(fill = treatment, y = Freq, x = subclass_main)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(values = c("#1e88e5", "#d81b60")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(breaks = c("cd4", "cd8", "n_k", "nkt",
                              "mon", "mac", "d_c", "mst", "b_c",
                              "tum", "ves"),
                   labels = paste0(c("cd4", "cd8", "n_k", "nkt",
                                     "mon", "mac", "d_c", "mst", "b_c",
                                     "tum", "ves"),
                                   " (n = ",
                                   c(9522, 8863, 2227, 3621,
                                     2721, 3380, 1485, 125, 783,
                                     1590, 309),
                                   ")")) +
  ylab("Fraction") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.height = unit(0, "cm"),
        legend.key.width = unit(0.25, "cm"))
dev.off()

# treatment group composition
tab <- data.frame(table(subclass_main = sce$subclass_main, treatment = sce$treatment))
tab$subclass_main <- factor(tab$subclass_main,
                               levels = c("cd4", "cd8", "n_k", "nkt",
                                          "mon", "mac", "d_c",
                                          "mst",
                                          "b_c",
                                          "tum", "ves"))

png("SFig1_tc.png", width = 600, height = 1450, res = 300)
ggplot(tab, aes(fill = subclass_main, y = Freq, x = treatment)) +
  geom_bar(position = "fill", stat = "identity", color = "black") +
  scale_fill_manual(breaks = c("cd4", "cd8", "n_k", "nkt",
                                "mon", "mac", "d_c",
                                "mst",
                                "b_c",
                                "tum", "ves"),
                     values = c("#ff9505", "#cc5803", "#e2711d", "#ffb627",
                                "#4cc9f0", "#4361ee", "#3a0ca3",
                                "#1b9aaa",
                                "#06d6a0",
                                "#424c55", "#d00000")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = c("control (n = 13422)",
                              "SBRT (n = 21204)")) +
  ylab("Fraction") +
  xlab("") +
  theme_classic() +
  theme(legend.position = "right",
        axis.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 0),
        legend.key.height = unit(0, "cm"),
        legend.key.width = unit(0.25, "cm"))
dev.off()

# subclassification TSNE
png(filename = "Fig2_TNSE_c.png", width = 1250, height = 1250, res = 300)
ggplot(a_c, aes(x = X1, y = X2,
                color = substr(subclass, 1, 3) %in% c("cd4", "cd8", "n_k", "nkt"))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey80", "darkgoldenrod2")) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none")
dev.off()

png(filename = "Fig5_TNSE_c.png", width = 1250, height = 1250, res = 300)
ggplot(a_c, aes(x = X1, y = X2,
                color = substr(subclass, 1, 3) %in% c("mon", "mac", "dc"))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("grey80", "dodgerblue3")) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none")
dev.off()

# track major genes from bulkRNA (Chow, PNAS 2020)
a_c <- cbind(a_c, t(logcounts(sce)[c("IFNG", "PDCD1", "IL16", "TNFSF10"), ]))

png("Fig1_bk.png", width = 600, height = 700, res = 300)
grid.arrange(
grobs = list(
ggplot(a_c, aes(x = subclass_main, y = IFNG)) +
  geom_jitter(alpha = 0.02, width = 0.2) +
  stat_summary(fun = "mean", geom = "point", color = "red", shape = 45, size = 10) +
  scale_x_discrete(limits = c("cd4", "cd8", "n_k", "nkt",
                              "mon", "mac", "d_c",
                              "mst", "b_c", "tum", "ves"),
                   labels = c("CD4", "CD8", "NK", "NKT",
                              "Mon", "Mac", "DC",
                              "Mast", "B cell", "Tum", "Ves")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,
                                    hjust = 0.5, face = "italic"),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"))
,
ggplot(a_c, aes(x = subclass_main, y = PDCD1)) +
  geom_jitter(alpha = 0.02, width = 0.2) +
  stat_summary(fun = "mean", geom = "point", color = "red", shape = 45, size = 10) +
  scale_x_discrete(limits = c("cd4", "cd8", "n_k", "nkt",
                              "mon", "mac", "d_c",
                              "mst", "b_c", "tum", "ves"),
                   labels = c("CD4", "CD8", "NK", "NKT",
                              "Mon", "Mac", "DC",
                              "Mast", "B cell", "Tum", "Ves")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,
                                    hjust = 0.5,
                                    face = "italic"),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"))
)
,
ncol = 1)

dev.off()

png("SFig1_bk.png", width = 600, height = 350, res = 300)
ggplot(a_c, aes(x = subclass_main, y = IL16)) +
  geom_jitter(alpha = 0.02, width = 0.2) +
  stat_summary(fun = "mean", geom = "point", color = "red", shape = 45, size = 10) +
  scale_x_discrete(limits = c("cd4", "cd8", "n_k", "nkt",
                              "mon", "mac", "d_c",
                              "mst", "b_c", "tum", "ves"),
                   labels = c("CD4", "CD8", "NK", "NKT",
                              "Mon", "Mac", "DC",
                              "Mast", "B cell", "Tum", "Ves")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,
                                    hjust = 0.5, face = "italic"),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"))
dev.off()

png("SFig_bk.png", width = 600, height = 350, res = 300)
ggplot(a_c, aes(x = subclass_main, y = TNFSF10)) +
  geom_jitter(alpha = 0.02, width = 0.2) +
  stat_summary(fun = "mean", geom = "point", color = "red", shape = 45, size = 10) +
  scale_x_discrete(limits = c("cd4", "cd8", "n_k", "nkt",
                              "mon", "mac", "d_c",
                              "mst", "b_c", "tum", "ves"),
                   labels = c("CD4", "CD8", "NK", "NKT",
                              "Mon", "Mac", "DC",
                              "Mast", "B cell", "Tum", "Ves")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7,
                                    hjust = 0.5, face = "italic"),
        axis.text.x = element_text(size = 7, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, color = "black"))
dev.off()
