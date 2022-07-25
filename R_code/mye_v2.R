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
library(limma)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)

# load data
mye <- readRDS("~/_Projects/Chow Nat Comm/data/mye_exp.RDS")

# load human gene labeling database
hs <- org.Hs.eg.db

# generate plot data.frame
pt <- cbind(x = reducedDim(mye, "TSNE")[, 1],
            y = reducedDim(mye, "TSNE")[, 2],
            data.frame(t(logcounts(mye)[c(
              "S100A8", "S100A9", "HLA-DRA",
              "CLEC10A", "CLEC9A", "FCER1A",
              "CD14", "FCGR3A", "CD68",
              "TNF", "IL1B", "CXCL8",
              "VEGFA", "MRC1", "FOLR2"
            ), ])),
            data.frame(colData(mye)))

#~~~MYE cluster subclassification~~~#
# cell subclassification
## assign subclass after examining variant and hallmark gene expression
pt$subclass <- ifelse(pt$hc_16 == 2, "d_c_1",
                      ifelse(pt$hc_16 == 3, "d_c_2",
                      ifelse(pt$hc_16 == 5, "d_c_3",
                      ifelse(pt$hc_16 == 10, "mon_1",
                      ifelse(pt$hc_16 == 4, "mon_2",
                      ifelse(pt$hc_16 == 1, "mon_3",
                      ifelse(pt$hc_16 == 11, "mon_4",
                      ifelse(pt$hc_16 == 6, "mon_5",
                      ifelse(pt$hc_16 == 8, "mon_6",
                      ifelse(pt$hc_16 == 12, "mac_1",
                      ifelse(pt$hc_16 == 16, "mac_2",
                      ifelse(pt$hc_16 == 9, "mac_3",
                      ifelse(pt$hc_16 == 13, "mac_4",
                      ifelse(pt$hc_16 == 7, "mac_5",
                      ifelse(pt$hc_16 == 15, "mac_6",
                      ifelse(pt$hc_16 == 14, "mac_7", NA
                             ))))))))))))))))

sub_nom <- data.frame(subclass = c(paste0("d_c_", 1:3),
                                   paste0("mac_", 1:7),
                                   paste0("mon_", 1:6)),
                      subclass_fine = c("DC_CD1C", "DC_C3", "DC_CLEC9A",
                                        "Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                                        "Mp_CD163", "Mp_SELENOP",
                                        "CD14_TNF", "CD14_IRF5",
                                        "CD16_CLEC7A", "CD16_IFITM3",
                                        "CD14_CD16_IRF3", "CD14_CD16_IL1B"))

pt <- left_join(pt, sub_nom, by = "subclass")

mye$subclass <- pt$subclass
mye$subclass_fine <- pt$subclass_fine

# saveRDS(mye, "~/_Projects/Chow Nat Comm/data/mye_FZ02.RDS")
# mye <- readRDS("~/_Projects/Chow Nat Comm/data/mye_FZ02.RDS")

sce <- SingleCellExperiment()
for(file_ in list.files("~/_Projects/Chow Nat Comm/data/sce_final", full.names = T)){
  if(nrow(sce) == 0){
    sce <- readRDS(file_)
  } else {
    sce <- cbind(sce, readRDS(file_))
  }
}

mye <- sce[, sce$class == "mye"]
reducedDim(mye, "UMAP") <- readRDS("~/_Projects/Chow Nat Comm/data/mye_UMAP.RDS")

sc_col <- colorRampPalette(c("#e4f9f5",
                             "#11999e",
                             "#005792",
                             "#00bbf0",
                             "#d9faff"))(16)

## subclustering of mye
pt <- cbind(x = reducedDim(mye, "UMAP")[, 1],
            y = reducedDim(mye, "UMAP")[, 2],
            data.frame(t(logcounts(mye)[c("HLA-DRA", "S100A8", "CD68", "FCER1A"), ])),
            data.frame(colData(mye)))

png(filename = "SFig6_UMAP_sc.png", width = 1100, height = 2000, res = 300)
grid.arrange(grobs = list(
  ggplot(pt, aes(x, y, color = subclass)) +
    scale_color_manual(values = sc_col, labels = 1:16) +
    geom_point(size = 0.1) +
    guides(color = guide_legend(override.aes = list(size = 1),
                                nrow = 2,
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
          legend.position = c(0.1, 0.1),
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

png(filename = "SFig6_sm.png", width = 1000, height = 1000, res = 300)
grid.arrange(grobs = lapply(
  c("HLA.DRA", "S100A8", "CD68", "FCER1A"),
  ts_plot)
  ,
  nrow = 2)
dev.off()

## plot on parent TSNE
########## by subclass ##########
pt <- cbind(x = reducedDim(mye, "TSNE")[, 1],
            y = reducedDim(mye, "TSNE")[, 2],
            data.frame(t(logcounts(mye)[c("S100A8", "S100A9", "HLA-DRA",
                                          "CLEC10A", "CLEC9A", "FCER1A",
                                          "CD14", "FCGR3A", "CD68",
                                          "TNF", "IL1B", "CXCL8",
                                          "VEGFA", "MRC1", "FOLR2"), ])),
            data.frame(colData(mye)))

labels <- data.frame(label = c("DC_CD1C", "DC_C3", "DC_CLEC9A",
                               "Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                               "Mp_CD163", "Mp_SELENOP",
                               "CD14_TNF", "CD14_IRF5",
                               "CD16_CLEC7A", "CD16_IFITM3",
                               "CD14_CD16_IRF3", "CD14_CD16_IL1B"),
                     x = c(11.5, 4, 11.5,
                           17, 12, 16, 5, -3,
                           -2, 0.5,
                           28, 25,
                           32, 35,
                           19.5, 16),
                     y = c(-16, -16, -12.25,
                           -24, -30, -18, -33.5, -22,
                           -29, -36.75,
                           -17, -21,
                           -4, -10,
                           -23, -33.5))

png(filename = "Fig5_TNSE_sc.png", width = 1000, height = 1000, res = 300)

ggplot(pt, aes(x, y)) +
  geom_point(size = 0.1,
             aes(color = subclass)) +
  scale_color_manual(values = sc_col) +
  scale_x_continuous(limits = c(-12, 37)) +
  scale_y_continuous(limits = c(-44, 0)) +
  geom_text_repel(data = filter(labels, grepl("SLC", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -1,
                  nudge_y = 5,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("C3|CLEC7A", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -7,
                  nudge_y = 2,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("CLEC9A", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 5,
                  nudge_y = 5) +
  geom_text_repel(data = filter(labels, grepl("CD1C", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -10,
                  nudge_y = 5) +
  geom_text_repel(data = filter(labels, grepl("IFITM3|IL1B", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 5,
                  nudge_y = -5) +
  geom_text_repel(data = filter(labels, grepl("IRF3", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 10,
                  nudge_y = -6) +
  geom_text_repel(data = filter(labels, grepl("IRF5|TNF", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 8,
                  nudge_y = -5) +
  geom_text_repel(data = filter(labels, grepl("CD80", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 3,
                  nudge_y = -7) +
  geom_text_repel(data = filter(labels, grepl("ACP", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 1,
                  nudge_y = 7,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("SELE", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -5,
                  nudge_y = -2,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("CD163", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = -4,
                  nudge_y = -3,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("IFNG", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 10,
                  nudge_y = -10,
                  min.segment.length = 0) +
  geom_text_repel(data = filter(labels, grepl("VEGFA", label)),
                  aes(label = label),
                  size = 2.5,
                  nudge_x = 0,
                  nudge_y = -6,
                  min.segment.length = 0) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "none")

dev.off()

########## by treatment ##########
png(filename = "Fig5_TSNE_tg.png", width = 1000, height = 1000, res = 300)

ggplot(pt, aes(x, y, color = treatment)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  guides(color = guide_legend(title = "Treatment",
                              title.position = "left",
                              nrow = 1,
                              override.aes = list(size = 2.5,
                                                  shape = 15))) +
  scale_x_continuous(limits = c(-12, 37)) +
  scale_y_continuous(limits = c(-44, 0)) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = c(0.35, 0.9),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))

dev.off()

########## plot of subclass by treatment ##########
png(filename = "Fig5_sc_tg.png", width = 500, height = 1250, res = 300)

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
  scale_y_discrete(limits = c(paste0("d_c_", 3:1), paste0("mac_", 7:1), paste0("mon_", 6:1)),
                   labels = c("DC_CLEC9A", "DC_C3", "DC_CD1C",  
                              "Mp_SELENOP", "Mp_CD163", "Mp_SLC40A1", "Mp_CD80", "Mp_ACP5",
                              "Mp_VEGFA", "Mp_IFNG", 
                              "CD14_CD16_IL1B", "CD14_CD16_IRF3", 
                              "CD16_IFITM3", "CD16_CLEC7A", 
                              "CD14_IRF5", "CD14_TNF")) +
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
ts_plot <- function(gene){
  ggplot(pt, aes(x, y)) +
    geom_point(data = pt[pt[, gene] == 0, ],
               color = "grey90",
               size = 0.1) +
    geom_point(data = pt[pt[, gene] > 0, ],
               aes_string(color = gene),
               size = 0.1) +
    scale_x_continuous(limits = c(-12, 37)) +
    scale_y_continuous(limits = c(-44, 0)) +
    scale_color_gradient(low = "grey90", high = "#003831") +
    ggtitle(gene) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "italic",
                                    size = 7,
                                    hjust = 0.5))
}

png(filename = "Fig5_sm.png", width = 1000, height = 1600, res = 300)

grid.arrange(grobs = lapply(
  c("S100A8", "S100A9", "HLA.DRA",
    "CLEC10A", "CLEC9A", "FCER1A",
    "CD14", "FCGR3A", "CD68",
    "TNF", "IL1B", "CXCL8",
    "VEGFA", "MRC1", "FOLR2"),
  ts_plot)
  ,
  ncol = 3)

dev.off()

## heatmap
### gene sets
lin <- c("CLEC10A", "FCER1A", "CD1C", "CLEC9A", # dc
         "S100A8", "S100A9", "SELL", "CD52", "FCN1",
         "CD14",
         "FCGR3A", "CDKN1C", "ITGAL", # mon
         "APOE", "APOC1", "C1QA", "C1QB", "C1QC") # mac
hla <- c(rownames(mye)[grepl("HLA-D", rownames(mye))])
d_c <- c("THBD", "THBS1", "C3")
mac <- c("CD80", "CD68", "CD86", "SPP1", "PDK4",
         "ACP5", "SLC40A1")
m_1 <- c("FCGR1A", "FCGR2A", "CXCL9", "CXCL10", "CXCL11",
         "MARCO", "CD40", "IDO1", "TNF", "IL1B",
         "IFNG", "IL6", "CXCL8", "SOCS3") # M1
m_2 <- c("VEGFA", "EREG", "FN1", "LGALS9", "CCL2", "MRC1",
         "CD163", "SELENOP", "TMEM176A", "TMEM176B",
         "FOLR2", "CLEC7A") # M2
ifn <- c("IRF8", "IRF3", "IRF5", "IFITM2", "IFITM3") # interferon response

### subset genes
sub <- mye[c(lin, hla, d_c, mac, m_1, m_2, ifn), ]

# clean up
rm(lin, hla, d_c, mac, m_1, m_2, ifn)

### calculate mean logcounts for genes of interest by cluster
heat_counts <- c()
for(x in sort(unique(mye$subclass))){
  a <- logcounts(sub)[, which(sub$subclass == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- c("DC_CD1C", "DC_C3", "DC_CLEC9A",
                           "Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                           "Mp_CD163", "Mp_SELENOP",
                           "CD14_TNF", "CD14_IRF5",
                           "CD16_CLEC7A", "CD16_IFITM3",
                           "CD14_CD16_IRF3", "CD14_CD16_IL1B")

### row annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             col = list(subclass = c("DC_CD1C" = sc_col[1],
                                                     "DC_C3" = sc_col[2],
                                                     "DC_CLEC9A" = sc_col[3],
                                                     "Mp_IFNG" = sc_col[4],
                                                     "Mp_VEGFA" = sc_col[5],
                                                     "Mp_ACP5" = sc_col[6],
                                                     "Mp_CD80" = sc_col[7],
                                                     "Mp_SLC40A1" = sc_col[8],
                                                     "Mp_CD163" = sc_col[9],
                                                     "Mp_SELENOP" = sc_col[10],
                                                     "CD14_TNF" = sc_col[11],
                                                     "CD14_IRF5" = sc_col[12],
                                                     "CD16_CLEC7A" = sc_col[13],
                                                     "CD16_IFITM3" = sc_col[14],
                                                     "CD14_CD16_IRF3" = sc_col[15],
                                                     "CD14_CD16_IL1B" = sc_col[16]
                             )),
                             which = "column",
                             height = unit(4, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

### col annotation relative abundance
set.seed(415); cont <- colnames(mye)[which(mye$treatment == "control")][sample(1:2686, 1500)]
set.seed(415); sbrt <- colnames(mye)[which(mye$treatment == "SBRT")][sample(1:4900, 1500)]

mye_rel <- mye[, c(cont, sbrt)]
mye_rel <- round(prop.table(table(mye_rel$subclass, mye_rel$treatment), 1) * 100, 1) %>%
  data.frame() %>%
  spread(Var2, Freq)
mye_rel <- HeatmapAnnotation("subclass\ncomposition" = anno_barplot(cbind(mye_rel$control, mye_rel$SBRT),
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
png(filename = "Fig5_hm.png", width = 1000, height = 3200, res = 300)

set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       top_annotation = col_ann,
                       bottom_annotation = mye_rel,
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
                       column_split = factor(c(rep("dc", 3),
                                               rep("mac", 7),
                                               rep("mon", 6)),
                                             levels = c("mon", "mac", "dc")),
                       row_split = factor(c(rep("lineage", 18),
                                            rep("HLA", 14),
                                            rep("DC", 3),
                                            rep("Mp", 7),
                                            rep("pro-inflammatory", 14),
                                            rep("anti-inflammatory", 12),
                                            rep("IFNG response", 5)),
                                          levels = c("lineage",
                                                     "HLA",
                                                     "IFNG response",
                                                     "DC",
                                                     "Mp",
                                                     "pro-inflammatory",
                                                     "anti-inflammatory")),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))

dev.off()

h_lgd = Legend(col_fun = col_fun, title = "Z score", at = c(-2, -1, 0, 1, 2),
               labels = c("< -2", "-1", "0", "1", "> 2"),
               direction = "horizontal",
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               legend_width = unit(3, "cm"))
b_lgd = Legend(labels = c("control", "SBRT"),
               title = "Treatment",
               legend_gp = gpar(fill = c("#1e88e5", "#d81b60")),
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               nrow = 1)
lgd <- packLegend(h_lgd, b_lgd, direction = "vertical")

png(filename = "Fig2_hm_lgd.png", width = 750, height = 250, res = 300)

draw(lgd)

dev.off()

### calculate mean logcounts for genes of interest by cluster separated by treatment group
sub$tc <- paste0(sub$treatment, " ", sub$subclass, " ", "hc:", sub$hc_16)

heat_counts <- c()
for(x in sort(unique(sub$tc))[c(1:7, 9:31)]){
  a <- logcounts(sub)[, which(sub$tc == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}
heat_counts <- cbind(heat_counts,
                     "control mac_6 hc:15" = logcounts(sub)[, which(sub$tc == "control mac_6 hc:15")])

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

heat_counts <- cbind(heat_counts,
                     "control mac_2 hc:16" = NA)

heat_counts <- heat_counts[, c(9:14, 4, 32, 5:7, 31, 8, 1:3, 25:30, 18:24, 15:17)]

colnames(heat_counts) <- paste(
  c(rep("control", 16), rep("SBRT", 16)),
  rep(c("CD14_TNF hc:10", "CD14_IRF5 hc:4",
        "CD16_CLEC7A hc:1", "CD16_IFITM3 hc:11",
        "CD14_CD16_IRF3 hc:6", "CD14_CD16_IL1B hc:8",
        "Mp_IFNG hc:12", "Mp_VEGFA hc:16", "Mp_ACP5 hc:9", "Mp_CD80 hc:13", "Mp_SLC40A1 hc:7",
        "Mp_CD163 hc:15", "Mp_SELENOP hc:14",
        "DC_CD1C hc:2", "DC_C3 hc:3", "DC_CLEC9A hc:5"), 2))

### row annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             treatment = c(rep("control", 16), rep("SBRT",16)),
                             col = list(subclass = c(
                               "control CD14_TNF hc:10" = sc_col[11],
                               "control CD14_IRF5 hc:4" = sc_col[12],
                               "control CD16_CLEC7A hc:1" = sc_col[13],
                               "control CD16_IFITM3 hc:11" = sc_col[14],
                               "control CD14_CD16_IRF3 hc:6" = sc_col[15],
                               "control CD14_CD16_IL1B hc:8" = sc_col[16],
                               "control Mp_IFNG hc:12" = sc_col[4],
                               "control Mp_VEGFA hc:16" = sc_col[5],
                               "control Mp_ACP5 hc:9" = sc_col[6],
                               "control Mp_CD80 hc:13" = sc_col[7],
                               "control Mp_SLC40A1 hc:7" = sc_col[8],
                               "control Mp_CD163 hc:15" = sc_col[9],
                               "control Mp_SELENOP hc:14" = sc_col[10],
                               "control DC_CD1C hc:2" = sc_col[1],
                               "control DC_C3 hc:3" = sc_col[2],
                               "control DC_CLEC9A hc:5" = sc_col[3],
                               "SBRT CD14_TNF hc:10" = sc_col[11],
                               "SBRT CD14_IRF5 hc:4" = sc_col[12],
                               "SBRT CD16_CLEC7A hc:1" = sc_col[13],
                               "SBRT CD16_IFITM3 hc:11" = sc_col[14],
                               "SBRT CD14_CD16_IRF3 hc:6" = sc_col[15],
                               "SBRT CD14_CD16_IL1B hc:8" = sc_col[16],
                               "SBRT Mp_IFNG hc:12" = sc_col[4],
                               "SBRT Mp_VEGFA hc:16" = sc_col[5],
                               "SBRT Mp_ACP5 hc:9" = sc_col[6],
                               "SBRT Mp_CD80 hc:13" = sc_col[7],
                               "SBRT Mp_SLC40A1 hc:7" = sc_col[8],
                               "SBRT Mp_CD163 hc:15" = sc_col[9],
                               "SBRT Mp_SELENOP hc:14" = sc_col[10],
                               "SBRT DC_CD1C hc:2" = sc_col[1],
                               "SBRT DC_C3 hc:3" = sc_col[2],
                               "SBRT DC_CLEC9A hc:5" = sc_col[3]
                               ),
                                        treatment = c("control" = "#1e88e5",
                                                      "SBRT" = "#d81b60")),
                             which = "column",
                             height = unit(3, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

########## plot heatmap ##########
png(filename = "SFig6_hm.png", width = 1450, height = 3250, res = 300)
set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       na_col = "grey80",
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
                       column_split = factor(c(rep("control", 16),
                                               rep("SBRT", 16)),
                                             levels = c("control", "SBRT")),
                       row_split = factor(c(rep("lineage", 18),
                                            rep("HLA", 14),
                                            rep("DC", 3),
                                            rep("Mp", 7),
                                            rep("pro-inflammatory", 14),
                                            rep("anti-inflammatory", 12),
                                            rep("IFNG response", 5)),
                                          levels = c("lineage",
                                                     "HLA",
                                                     "IFNG response",
                                                     "DC",
                                                     "Mp",
                                                     "pro-inflammatory",
                                                     "anti-inflammatory")),
                       row_gap = unit(2, "mm"),
                       column_gap = unit(5, "mm"),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))
dev.off()

## clean up
rm(list = ls()[!ls() %in% c("mye", "hs", "sc_col", "col_ann")])

##########~~~Monocyte/Macrophage analysis~~~##########

## dge by limma
mnm <- mye[, grepl("m", mye$subclass)]
mnm <- mnm[rowSums(counts(mnm)) > 0, ] # only genes detectable within this subset
counts <- log2(1 + calculateCPM(mnm))
design <- model.matrix(~I(treatment == "SBRT"), colData(mnm))
l <- lmFit(counts, design)
e <- eBayes(l)
deg <- topTable(e, n = Inf)
deg$sig <- ifelse(abs(deg$logFC) < 1 |
                    deg$adj.P.Val > 0.05,
                  "ns",
                  ifelse(deg$logFC > 1,
                         "SBRT",
                         "cont."))

rm(counts, design, l, e)

########## MnM volcano plot ##########
deg_ <- deg
deg_$adj.P.Val[deg_$adj.P.Val < 10^-300] <- 10^-299.5 # rescale for cosmetics

png(filename = "Fig6_vp.png", width = 1000, height = 850, res = 300)

ggplot(deg_, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.5, aes(color = sig)) +
  scale_color_manual(breaks = c("cont.", "ns", "SBRT"),
                     values = c("#1e88e5", "grey", "#d81b60")) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  shape = 15),
                              nrow = 1)) +
  labs(color = "DGE") +
  xlab(expression("control    " %<-% "    log2FC    " %->% "    SBRT")) +
  geom_text_repel(data = filter(deg_, ID %in% c("VEGFA")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -3,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg_, ID %in% c("HIF1A")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 1,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg_, ID %in% c("THBS1", "MIF")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 1.5,
                  nudge_y = -8) +
  geom_text_repel(data = filter(deg_, ID %in% c("CXCL8")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 1.5,
                  nudge_y = -20) +
  geom_text_repel(data = filter(deg_, ID %in% c("FCGR3A")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 2,
                  nudge_y = -14) +
  geom_text_repel(data = filter(deg_, ID %in% c("IRF8", "FOLR2")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 4,
                  nudge_y = 2) +
  geom_text_repel(data = filter(deg_, ID %in% c("THBD", "IL1B", "CD14")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 1.5,
                  nudge_y = 5) +
  geom_text_repel(data = filter(deg_, ID %in% c("APOE")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 1.5,
                  nudge_y = 15) +
  geom_text_repel(data = filter(deg_, ID %in% c("HLA-DRA")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -1.5,
                  nudge_y = 40) +
  geom_text_repel(data = filter(deg_, ID %in% c("FCN1")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -1,
                  nudge_y = 20) +
  geom_text_repel(data = filter(deg_, ID %in% c("CD52")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -2,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg_, ID %in% c("CDKN1C")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -1,
                  nudge_y = -1) +
  scale_x_continuous(limits = c(-5.6, 5.6),
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

### probe lists
m_1 <- c("FCGR1A", "FCGR2A", "CXCL9", "CXCL10", "CXCL11",
         "MARCO", "CD40", "IDO1", "TNF", "IL1B",
         "IFNG", "IL6", "CXCL8", "SOCS3", "IRF5",
         "CD80", "CD86", "CD40", "IRF1", "IRF5", "IRF8") # M1
m_1 <- AnnotationDbi::select(hs, keys = m_1, columns = c("ENTREZID"), keytype = "SYMBOL") # match gene symbol to ENTREZID
m_1 <- m_1[which(is.na(m_1$ENTREZID) == F), "ENTREZID"] # filter genes with missing ENTREZID

m_2 <- c("VEGFA", "EREG", "FN1", "CCL2", "MRC1",
         "CD163", "SELENOP", "TMEM176A", "TMEM176B",
         "FOLR2", "CLEC7A", "MSR1", "VEGFB", "CD68",
         "MMP14", "MMP9", "MMP19") # M2
m_2 <- AnnotationDbi::select(hs, keys = m_2, columns = c("ENTREZID"), keytype = "SYMBOL") # match gene symbol to ENTREZID
m_2 <- m_2[which(is.na(m_2$ENTREZID) == F), "ENTREZID"] # filter genes with missing ENTREZID

mac <- c("C1QA", "C1QB", "C1QC",
         "APOE", "APOC1", "HLA-DRA") # mac
mac <- AnnotationDbi::select(hs, keys = mac, columns = c("ENTREZID"), keytype = "SYMBOL") # match gene symbol to ENTREZID
mac <- mac[which(is.na(mac$ENTREZID) == F), "ENTREZID"] # filter genes with missing ENTREZID

### calculate GSEA of lists
score_lists <- list("Pro-inflammatory" = m_1,
                    "Anti-inflammatory" = m_2,
                    "Macrophage" = mac) # rename "Pathway" to name of choice

fgsea(pathways = score_lists,
      stats    = ent_ord,
      minSize  = 1,
      maxSize  = 500)

statsAdj <- ent_ord / max(abs(ent_ord))

########## Pro-inflammatory GSEA ##########
pathway <- unname(as.vector(na.omit(match(score_lists[[1]], names(statsAdj))))) # matches query list to ordered list
pathway <- sort(pathway) # orders query list to ordered list

gseaRes <- calcGseaStat(statsAdj,
                        selectedStats = pathway,
                        returnAllExtremes = TRUE) # calculate GSEA stats

bottoms <- gseaRes$bottoms # extract bottom for RES
tops <- gseaRes$tops # extract top for RES

n <- length(statsAdj) # to make ggplot data.frame
xs <- as.vector(rbind(pathway - 1, pathway)) # x coordinate
ys <- as.vector(rbind(bottoms, tops)) # y coordinate

toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0)) # add additional coordinates to end at y = 0

diff <- (max(tops) - min(bottoms)) / 8 # scales query line height

png(filename = "Fig6_pi.png", width = 700, height = 300, res = 300)

ggplot(toPlot, aes(x = x, y = y)) +
  xlab("Gene rank") +
  ylab("Enrichment score") +
  ggtitle("Pro-inflammatory, NES = 1.33, padj < 0.05") +
  geom_point(color = "green", size = 0.1) +
  geom_hline(yintercept = max(tops), color = "red", linetype = "dashed") +
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") + theme_classic() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff,
                             xend = x, yend = diff),
               size = 0.25) +
  theme(axis.text.x = element_text(size = 7,
                                   color = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.title.x = element_text(size = 7,
                                    color = "white"),
        axis.title.y = element_text(size = 7,
                                    color = "white"),
        plot.title = element_text(size = 7))

dev.off()

########## Anti-inflammatory GSEA ##########
pathway <- unname(as.vector(na.omit(match(score_lists[[2]], names(statsAdj))))) # matches query list to ordered list
pathway <- sort(pathway) # orders query list to ordered list

gseaRes <- calcGseaStat(statsAdj,
                        selectedStats = pathway,
                        returnAllExtremes = TRUE) # calculate GSEA stats

bottoms <- gseaRes$bottoms # extract bottom for RES
tops <- gseaRes$tops # extract top for RES

n <- length(statsAdj) # to make ggplot data.frame
xs <- as.vector(rbind(pathway - 1, pathway)) # x coordinate
ys <- as.vector(rbind(bottoms, tops)) # y coordinate

toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0)) # add additional coordinates to end at y = 0

diff <- (max(tops) - min(bottoms)) / 8 # scales query line height

png(filename = "Fig6_ai.png", width = 700, height = 300, res = 300)

ggplot(toPlot, aes(x = x, y = y)) +
  xlab("Gene rank") +
  ylab("Enrichment score") +
  ggtitle("Anti-inflammatory, NES = 1.37, padj < 0.05") +
  geom_point(color = "green", size = 0.1) +
  geom_hline(yintercept = max(tops), color = "red", linetype = "dashed") +
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") + theme_classic() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff,
                             xend = x, yend = diff),
               size = 0.25) +
  theme(axis.text.x = element_text(size = 7,
                                   color = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.title.x = element_text(size = 7,
                                    color = "white"),
        axis.title.y = element_text(size = 7,
                                    color = "black"),
        plot.title = element_text(size = 7))

dev.off()

########## Macrophage GSEA ##########
pathway <- unname(as.vector(na.omit(match(score_lists[[3]], names(statsAdj))))) # matches query list to ordered list
pathway <- sort(pathway) # orders query list to ordered list

gseaRes <- calcGseaStat(statsAdj,
                        selectedStats = pathway,
                        returnAllExtremes = TRUE) # calculate GSEA stats

bottoms <- gseaRes$bottoms # extract bottom for RES
tops <- gseaRes$tops # extract top for RES

n <- length(statsAdj) # to make ggplot data.frame
xs <- as.vector(rbind(pathway - 1, pathway)) # x coordinate
ys <- as.vector(rbind(bottoms, tops)) # y coordinate

toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0)) # add additional coordinates to end at y = 0

diff <- (max(tops) - min(bottoms)) / 8 # scales query line height

png(filename = "Fig6_mp.png", width = 700, height = 300, res = 300)

ggplot(toPlot, aes(x = x, y = y)) +
  xlab(expression("SBRT    " %<-% "    Gene rank    " %->% "    control")) +
  ylab("Enrichment score") +
  ggtitle("Macrophage, NES = 1.33, padj < 0.05") +
  scale_y_continuous(breaks = seq(0, 0.5, 0.5)) +
  geom_point(color = "green", size = 0.1) +
  geom_hline(yintercept = max(tops), color = "red", linetype = "dashed") +
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") + theme_classic() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff,
                             xend = x, yend = diff),
               size = 0.25) +
  theme(axis.text.x = element_text(size = 7,
                                   color = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.title.x = element_text(size = 7,
                                    color = "black"),
        axis.title.y = element_text(size = 7,
                                    color = "white"),
        plot.title = element_text(size = 7))

dev.off()

## pseudotime
### recalculate UMAP for Mo and Mp only
set.seed(415); mnm <- runPCA(mnm, ntop = 125)
set.seed(415); mnm <- runUMAP(mnm, dimred = "PCA", n_dimred = c(1:5, 7:50)) # remove PCs to plot clusters together
mnm$hc <- factor(cutree(hclust(dist(reducedDim(mnm, "UMAP"))), 6)) # remove stray cluster
mnm <- mnm[, !mnm$hc == 6]

### calculate UMAP, use existing subcluster id
mnm <- slingshot(mnm, clusterLabels = "hc", reducedDim = "UMAP", start.clus = "2")

mnm$Pseudotime <- ifelse(is.na(mnm$slingPseudotime_1),
                         mnm$slingPseudotime_2,
                         mnm$slingPseudotime_1)

pt <- cbind(x = reducedDims(mnm)$UMAP[, 1],
            y = -reducedDims(mnm)$UMAP[, 2],
            data.frame(colData(mnm)[, c("treatment", "Pseudotime", "slingPseudotime_1", "slingPseudotime_2", "subclass", "hc")]),
            data.frame(t(logcounts(mnm)[c(
              "CD52", "S100A8", "S100A9",
              "CD14", "FCGR3A", "CD68", "CD80", "CD86", "CD40", "HLA-DRA",
              "C1QC", "APOE", "APOC1",
              "IL1B", "CXCL8", "FOLR2", "CD163", "VEGFA", "TNF", "SELENOP",
              "IRF8", "IRF3", "IRF5", "IFITM2", "IFITM3", "IRF1"
            ), ])))

### calculate lineage
lin <- getLineages(reducedDims(mnm)$UMAP, mnm$hc, start.clus = "2")
cur <- getCurves(lin)

cur1 <- data.frame(x = cur@metadata$curves$Lineage1$s[, 1],
                   y = -cur@metadata$curves$Lineage1$s[, 2])
cur2 <- data.frame(x = cur@metadata$curves$Lineage2$s[, 1],
                   y = -cur@metadata$curves$Lineage2$s[, 2])

# clean up
rm(lin, cur)

########## MnM pseudotime ##########
png(filename = "Fig6_ps.png", width = 1450, height = 900, res = 300)

grid.arrange(grobs = list(
  ggplot(pt, aes(x, y)) +
    geom_point(aes(color = slingPseudotime_2),
               size = 0.1) +
    scale_color_distiller(palette = "BuGn", direction = 1, limits = c(0, 8), na.value = "grey75") +
    guides(color = guide_colorbar(barwidth = 5, barheight = 0.25,
                                  frame.colour = "black", ticks.colour = "black",
                                  title = "Non-classical Mo\nlineage pseudotime",
                                  title.position = "left")) +
    geom_smooth(data = cur1,
                size = 0.5,
                color = "black",
                se = F) +
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
    scale_color_distiller(palette = "BuPu", direction = 1, limits = c(0, 16), na.value = "grey75") +
    guides(color = guide_colorbar(barwidth = 5, barheight = 0.25,
                                  frame.colour = "black", ticks.colour = "black",
                                  title = "Macrophage lineage\npseudotime",
                                  title.position = "left")) +
    geom_smooth(data = cur1,
                size = 0.5,
                color = "black",
                se = F) +
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
),
layout_matrix = cbind(c(1, 2)))

dev.off()

########## MnM by treatment group ##########
png(filename = "Fig6_tg.png", width = 1450, height = 800, res = 300)

grid.arrange(grobs = list(
  ggplot(pt, aes(x, y)) +
    ggtitle("control\n") +
    geom_vline(xintercept = -5.4, color = "darkorchid4", linetype = "dotdash") +
    geom_point(data = filter(pt, treatment == "SBRT"),
               color = "grey75",
               size = 0.1) +
    geom_point(data = filter(pt, treatment == "control"),
               color = "#1e88e5",
               size = 0.1) +
    geom_smooth(data = cur1,
                size = 0.5,
                color = "black",
                se = F) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          plot.title = element_text(size = 7,
                                    hjust = 0.5))
  ,
  ggplot(pt, aes(x, y)) +
    ggtitle("SBRT\n") +
    geom_vline(xintercept = -5.4, color = "darkorchid4", linetype = "dotdash") +
    geom_point(data = filter(pt, treatment == "control"),
               color = "grey75",
               size = 0.1) +
    geom_point(data = filter(pt, treatment == "SBRT"),
               color = "#d81b60",
               size = 0.1) +
    geom_smooth(data = cur1,
                size = 0.5,
                color = "black",
                se = F) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          plot.title = element_text(size = 7,
                                    hjust = 0.5))
  ,
  ggplot(pt, aes(x, fill = treatment)) +
    geom_density() +
    geom_vline(xintercept = -5.4, color = "darkorchid4", linetype = "dotdash") +
    xlab("Pseudotime") +
    ylab("Cell density") +
    scale_x_continuous(expand = c(0, 0),) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.2)) +
    scale_fill_manual(values = c("#1e88e5", "#d81b60")) +
    theme_classic() +
    theme(aspect.ratio = 0.5,
          legend.position = "none",
          plot.title = element_text(size = 7,
                                    color = "white",
                                    hjust = 0.5),
          axis.title.x = element_text(size = 7,
                                      color = "black"),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 7,
                                      color = "black"),
          axis.text.y = element_text(size = 7,
                                     color = "black"),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 7,
                                    color = "white",
                                    hjust = 0.5),
          strip.background = element_blank()) +
    facet_wrap(~treatment, nrow = 2, scales = "free")
), layout_matrix = rbind(c(1, 1, 3, 3),
                         c(2, 2, 3, 3)))

dev.off()

########## MnM gene expression by UMAP ##########
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
    geom_smooth(data = cur1,
                size = 0.5,
                color = "black",
                se = F) +
    geom_point(data = cur2,
               size = 0.1,
               shape = 46) +
    theme_void() +
    theme(aspect.ratio = 0.33,
          legend.title = element_text(size = 7, hjust = 0, face = "italic"),
          legend.text = element_text(size = 7),
          legend.position = "bottom")
}

png(filename = "Fig6_ps_mnm.png", width = 1450, height = 750, res = 300)

grid.arrange(grobs = lapply(
  c("CD52", "S100A9", "HLA.DRA",
    "CD14", "FCGR3A", "CD68",
    "C1QC", "APOE", "APOC1"),
  ps_plot
), nrow = 3)

dev.off()

########## mon/mac chisq ##########
set.seed(415); cont <- colnames(mnm)[which(mnm$treatment == "control")][sample(1:1702, 1500)]
set.seed(415); sbrt <- colnames(mnm)[which(mnm$treatment == "SBRT")][sample(1:4245, 1500)]
mnm_1500 <- mnm[, c(cont, sbrt)]

chisq.test(table(mnm_1500$treatment, substr(mnm_1500$subclass, 1, 3)))
cs <- data.frame(table(treatment = mnm_1500$treatment, subclass = substr(mnm_1500$subclass, 1, 3)))
cs$x0 <- c(2, 2, 1, 1)
cs$y0 <- c(2, 1, 2, 1)
cs$frac <- cs$Freq/max(cs$Freq)
cs$Freq <- cs$Freq/sum(cs$Freq)

png(filename = "Fig6_ps_cs.png", width = 800, height = 480, res = 300)

ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = frac/2, fill = Freq), data = cs) +
  scale_fill_distiller(palette = 7, direction = 1) +
  guides(fill = guide_colorbar(title = "Fraction\nof cells",
                               title.position = "top",
                               barheight = 3.5,
                               barwidth = 0.5,
                               frame.colour = "black",
                               ticks.colour = "black")) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("monocyte", "macrophage")) +
  scale_y_continuous(limits = c(0.5, 2.5), expand = c(0,0), breaks = c(1, 2), labels = c("SBRT", "control")) +
  theme_minimal() +
  theme(aspect.ratio = 1,
        legend.position = "right",
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

########## macrophage pseudotime smoothed expression by treatment group ##########
png(filename = "Fig6_ps_sm.png", width = 2000, height = 500, res = 300)

pt_tg <- filter(pt, !is.na(slingPseudotime_1))

pt_tg <- gather(pt_tg,
                "IL1B", "CXCL8", "FOLR2", "SELENOP", "VEGFA", "TNF",
                "CD68", "HLA.DRA",
                key = "gene", value = "logcounts")
pt_tg$gene <- factor(pt_tg$gene,
                     levels = c("HLA.DRA", "IL1B", "CXCL8", "TNF",
                                "CD68", "VEGFA", "FOLR2", "SELENOP"))

ggplot(pt_tg,
       aes(y = logcounts)) +
  geom_smooth(aes(x = slingPseudotime_1, color = treatment), se = F) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Macophage lineage pseudotime") +
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
             scale = "free_y",
             nrow = 2)
dev.off()

########## macrophage responses ##########
# show Mp in pseudotime
pt$mye_type <- ifelse(pt$subclass %in% c(paste0("mac_", 1:7)), "Mac",
                      ifelse(pt$subclass %in% c("mon_1", "mon_2"), "CD14",
                             ifelse(pt$subclass %in% c("mon_3", "mon_4"), "CD16", "CD14_CD16")))

png(filename = "SFig7_ps_mp.png", width = 1000, height = 333, res = 300)
ggplot(pt, aes(x, y, color = mye_type)) +
  scale_color_manual(values = c("dodgerblue", "cyan", "blue2", "deepskyblue3")) +
  geom_point(size = 0.1) +
  geom_text_repel(data = data.frame(x = 5.5, y = 0.4),
                  label = "Macrophage",
                  color = "black",
                  nudge_x = 1.25,
                  nudge_y = 1,
                  size = 2.25,
                  min.segment.length = 0) +
  geom_text_repel(data = data.frame(x = 0, y = 1.5),
                  label = "Mon: CD14+ CD16+",
                  color = "black",
                  nudge_x = 0,
                  nudge_y = 1,
                  size = 2.25,
                  min.segment.length = 0) +
  geom_text_repel(data = data.frame(x = -5, y = -1),
                  label = "Mon: CD14+",
                  color = "black",
                  nudge_x = 1,
                  nudge_y = 1,
                  size = 2.25,
                  min.segment.length = 0) +
  geom_text_repel(data = data.frame(x = -6.5, y = 2.5),
                  label = "Mon: CD16+",
                  color = "black",
                  nudge_x = 1,
                  nudge_y = 1,
                  size = 2.25,
                  min.segment.length = 0) +
  theme_void() +
  theme(legend.position = "none",
        aspect.ratio = 0.33)
dev.off()

## dge by limma
mac <- mye[, grepl("mac", mye$subclass)]
mac <- mac[rowSums(counts(mac)) > 0, ] # only genes detectable within this subset
counts <- log2(1 + calculateCPM(mac))
design <- model.matrix(~I(treatment == "SBRT"), colData(mac))
l <- lmFit(counts, design)
e <- eBayes(l)
deg <- topTable(e, n = Inf)
deg$sig <- ifelse(abs(deg$logFC) < 1 |
                    deg$adj.P.Val > 0.05,
                  "ns",
                  ifelse(deg$logFC > 1,
                         "SBRT",
                         "cont."))
rm(counts, design, l, e)

### Interferon sensitive
deg_sub <- filter(deg, grepl("^IRF|^IFI", ID), logFC != 0)
deg_sub <- deg_sub[order(deg_sub$logFC, decreasing = T), ]

### subset genes
sub <- mac[deg_sub$ID, ]

### calculate mean logcounts for genes of interest by cluster
heat_counts <- c()
for(x in sort(unique(sub$subclass))){
  a <- logcounts(sub)[, which(sub$subclass == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- c("Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                           "Mp_CD163", "Mp_SELENOP")

### col annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             col = list(subclass = c("Mp_IFNG" = sc_col[4],
                                                     "Mp_VEGFA" = sc_col[5],
                                                     "Mp_ACP5" = sc_col[6],
                                                     "Mp_CD80" = sc_col[7],
                                                     "Mp_SLC40A1" = sc_col[8],
                                                     "Mp_CD163" = sc_col[9],
                                                     "Mp_SELENOP" = sc_col[10]
                             )),
                             which = "column",
                             height = unit(4, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

### row annotation relative abundance
row_deg <- HeatmapAnnotation("log2FC" = anno_barplot(deg_sub$logFC,
                                                    gp = gpar(fill = c(rep("#d81b60", 14),
                                                                       rep("grey80", 16))),
                                                    bar_width = 0.9),
                             which = "row",
                             width = unit(2, "cm"),
                             annotation_name_gp = gpar(fontsize = 7),
                             annotation_name_rot = 0)

### heatmap scale
col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                colorRampPalette(c("purple", "black", "yellow"))(5))

########## plot heatmap ##########
png(filename = "SFig7_ifrs.png", width = 700, height = 1500, res = 300)
set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       bottom_annotation = col_ann,
                       left_annotation = row_deg,
                       show_heatmap_legend = F,
                       cluster_rows = F,
                       cluster_columns = T,
                       show_row_dend = F,
                       show_row_names = T,
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title_gp = gpar(fontsize = 7),
                       show_column_dend = T,
                       show_column_names = T,
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 7),
                       column_title_gp = gpar(fontsize = 0),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))
dev.off()

### Hypoxia sensitive
deg_sub <- filter(deg, grepl("^HIF|^HIG|EPAS1|VEGF", ID), logFC != 0)
deg_sub <- deg_sub[order(deg_sub$logFC, decreasing = T), ]

### subset genes
sub <- mac[deg_sub$ID, ]

### calculate mean logcounts for genes of interest by cluster
heat_counts <- c()
for(x in sort(unique(sub$subclass))){
  a <- logcounts(sub)[, which(sub$subclass == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- c("Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                           "Mp_CD163", "Mp_SELENOP")

### row annotation relative abundance
row_deg <- HeatmapAnnotation("log2FC" = anno_barplot(deg_sub$logFC,
                                                    gp = gpar(fill = c(rep("#d81b60", 5),
                                                                       rep("grey80", 7))),
                                                    bar_width = 0.9),
                             which = "row",
                             width = unit(2, "cm"),
                             annotation_name_gp = gpar(fontsize = 7),
                             annotation_name_rot = 0)

########## plot heatmap ##########
png(filename = "SFig7_hpox.png", width = 700, height = 900, res = 300)
set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       bottom_annotation = col_ann,
                       left_annotation = row_deg,
                       show_heatmap_legend = F,
                       cluster_rows = F,
                       cluster_columns = T,
                       show_row_dend = F,
                       show_row_names = T,
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title_gp = gpar(fontsize = 7),
                       show_column_dend = T,
                       show_column_names = T,
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 7),
                       column_title_gp = gpar(fontsize = 0),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))
dev.off()

### Mac co-stim and co-inhibitor
deg_cs <- filter(deg, ID %in% c("CD80", "CD86", "CD70", "CD40",
                                "TNFSF14", "TNFSF9", "TNFSF4"), logFC != 0)
deg_cs <- deg_cs[order(deg_cs$logFC, decreasing = T), ]

deg_ci <- filter(deg, ID %in% c("TNFRSF14", "CD274", "LGALS9",
                                "CD48", "TEK", "CCL2", "IL6"), logFC != 0)
deg_ci <- deg_ci[order(deg_ci$logFC, decreasing = T), ]

deg_sub <- rbind(deg_cs, deg_ci)

### subset genes
sub <- mac[deg_sub$ID, ]

### calculate mean logcounts for genes of interest by cluster
heat_counts <- c()
for(x in sort(unique(sub$subclass))){
  a <- logcounts(sub)[, which(sub$subclass == x)]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  heat_counts <- cbind(heat_counts, a)
}

heat_counts <- t(scale(t(heat_counts)))
heat_counts[heat_counts > 2] <- 2
heat_counts[heat_counts < -2] <- -2

colnames(heat_counts) <- c("Mp_IFNG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1",
                           "Mp_CD163", "Mp_SELENOP")

### col annotation
col_ann <- HeatmapAnnotation(subclass = colnames(heat_counts),
                             col = list(subclass = c("Mp_IFNG" = sc_col[4],
                                                     "Mp_VEGFA" = sc_col[5],
                                                     "Mp_ACP5" = sc_col[6],
                                                     "Mp_CD80" = sc_col[7],
                                                     "Mp_SLC40A1" = sc_col[8],
                                                     "Mp_CD163" = sc_col[9],
                                                     "Mp_SELENOP" = sc_col[10]
                             )),
                             which = "column",
                             height = unit(4, "mm"),
                             gp = gpar(col = "white"),
                             show_legend = F,
                             annotation_name_gp = gpar(fontsize = 0))

### row annotation relative abundance
row_deg <- HeatmapAnnotation("log2FC" = anno_barplot(deg_sub$logFC,
                                                     gp = gpar(fill = c("#d81b60",
                                                                        rep("grey80", 6),
                                                                        "#d81b60",
                                                                        rep("grey80", 6))),
                                                     bar_width = 0.9),
                             which = "row",
                             width = unit(2, "cm"),
                             annotation_name_gp = gpar(fontsize = 7),
                             annotation_name_rot = 0)

########## plot heatmap ##########
png(filename = "SFig7_cos.png", width = 720, height = 1000, res = 300)
set.seed(415); Heatmap(heat_counts,
                       col = col_fun,
                       border = T,
                       bottom_annotation = col_ann,
                       left_annotation = row_deg,
                       show_heatmap_legend = F,
                       cluster_rows = F,
                       cluster_columns = T,
                       show_row_dend = F,
                       show_row_names = T,
                       row_names_side = "right",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title_gp = gpar(fontsize = 7),
                       row_split = factor(c(rep("Co-stimulator", 7),
                                            rep("Co-inhibitor", 7)),
                                          levels = c("Co-stimulator",
                                                     "Co-inhibitor")),
                       show_column_dend = T,
                       show_column_names = T,
                       column_names_side = "bottom",
                       column_names_gp = gpar(fontsize = 7),
                       column_title_gp = gpar(fontsize = 0),
                       width = ncol(heat_counts)*unit(3, "mm"), 
                       height = nrow(heat_counts)*unit(3, "mm"))
dev.off()

# legend
h_lgd = Legend(col_fun = col_fun, title = "Z score", at = c(-2, -1, 0, 1, 2),
               labels = c("< -2", "-1", "0", "1", "> 2"),
               direction = "vertical",
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               legend_height = unit(2, "cm"),
               grid_width = unit(0.25, "cm"))
g_lgd = Legend(labels = c("ns", "SBRT"),
               title = "DEG",
               legend_gp = gpar(fill = c("grey80", "#d81b60")),
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               grid_width = unit(0.25, "cm"),
               grid_height = unit(0.25, "cm"),
               nrow = 2)
lgd <- packLegend(h_lgd, g_lgd, direction = "vertical", column_gap = unit(5, "mm"))

png(filename = "SFig7_hm_lgd.png", width = 150, height = 450, res = 300)

draw(lgd)

dev.off()

h_lgd = Legend(col_fun = col_fun, title = "Z score", at = c(-2, -1, 0, 1, 2),
               labels = c("< -2", "-1", "0", "1", "> 2"),
               direction = "horizontal",
               title_gp = gpar(fontsize = 7),
               title_position = "leftcenter",
               labels_gp = gpar(fontsize = 7),
               legend_height = unit(2, "cm"),
               grid_height = unit(0.25, "cm"))
g_lgd = Legend(labels = c("ns", "SBRT"),
               title = "    DEG",
               legend_gp = gpar(fill = c("grey80", "#d81b60")),
               title_gp = gpar(fontsize = 7),
               title_position = "leftcenter",
               labels_gp = gpar(fontsize = 7),
               grid_width = unit(0.25, "cm"),
               grid_height = unit(0.25, "cm"),
               nrow = 1)
lgd <- packLegend(g_lgd, h_lgd, direction = "vertical", row_gap = unit(7, "mm"))

png(filename = "SFig7_hm_hlgd.png", width = 1000, height = 150, res = 300)

draw(lgd)

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
  filter(grepl("Interferon|ypox", Description))
sub_GSEA$Description <- factor(sub_GSEA$Description,
                               levels = sub_GSEA$Description[order(sub_GSEA$NES, decreasing = F)])

########## Mac GSEA ##########
png(filename = "SFig7_gsea.png", width = 910, height = 520, res = 300)

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
  ggtitle("Reactome Pathway Enrichment") +
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

