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
library(limma)
library(dplyr)
library(tidyr)
library(ReactomePA)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)

# load data
tum <- readRDS("~/_Projects/Chow Nat Comm/data/tum_FZ02.RDS")
all <- readRDS("~/_Projects/Chow Nat Comm/data/sce_classed_FZ02.RDS")

# load human gene labeling database
hs <- org.Hs.eg.db

# generate plot data.frame
pt <- cbind(x = reducedDim(all, "TSNE")[, 1],
            y = reducedDim(all, "TSNE")[, 2],
            data.frame(colData(all)))

## plot on parent TSNE
########## by treatment ##########
png(filename = "Fig7_TSNE.png", width = 600, height = 600, res = 300)

ggplot() +
  geom_point(data = filter(pt, class == "tum"), aes(x, y), color = "#424c55", size = 0.1) +
  geom_point(data = filter(pt, class != "tum"), aes(x, y), color = "grey75", size = 0.1) +
  geom_text_repel(aes(x = -25, y = -20, label = "tumor"),
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -3,
                  nudge_y = -8) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))

dev.off()

# clean up
rm(all)

pt <- cbind(x = reducedDim(tum, "TSNE")[, 1],
            y = reducedDim(tum, "TSNE")[, 2],
            data.frame(t(logcounts(tum)[c(
              "HLA-A", "HLA-B", "HLA-C",
              "IFNGR1", "CA9", "CALR",
              "FAS", "BAD", "BAX"
            ), ])),
            data.frame(colData(tum)))

########## by treatment ##########
png(filename = "Fig7_tum_TSNE.png", width = 750, height = 750, res = 300)

ggplot(pt, aes(x, y, color = treatment)) +
  geom_point(size = 0.1) +
  ylim(c(-29, -11)) +
  xlim(c(-33, -15)) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  guides(color = guide_legend(title = "Treatment",
                              title.position = "left",
                              nrow = 1,
                              override.aes = list(size = 2.5,
                                                  shape = 15))) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7))

dev.off()

########## small multiples ##########
# generate plot data.frame
ts_plot <- function(gene){
  ggplot(pt, aes(x, y)) +
    geom_point(data = pt[pt[, gene] == 0, ],
               color = "grey90",
               size = 0.1) +
    geom_point(data = pt[pt[, gene] > 0, ],
               aes_string(color = gene),
               size = 0.1) +
    ylim(c(-29, -11)) +
    xlim(c(-33, -15)) +
    scale_color_gradient(low = "grey90", high = "#003831") +
    ggtitle(gene) +
    theme_void() +
    theme(legend.position = "none",
          aspect.ratio = 1,
          plot.title = element_text(face = "italic",
                                    size = 7,
                                    hjust = 0.5))
}

png(filename = "Fig7_sm.png", width = 1000, height = 900, res = 300)

grid.arrange(grobs = lapply(
  c("HLA.A", "HLA.B", "HLA.C",
    "IFNGR1", "CA9", "CALR",
    "FAS", "BAD", "BAX"),
  ts_plot)
  ,
  ncol = 3)

dev.off()

##########~~~Tumor analysis~~~##########
## dge by limma
counts <- log2(1 + calculateCPM(tum))
tum <- tum[rowSums(counts(tum)) > 0, ] # only genes detectable within this subset
design <- model.matrix(~I(treatment == "SBRT"), colData(tum))
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

########## tumor volcano plot ##########
png(filename = "Fig7_vp.png", width = 700, height = 750, res = 300)

ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 0.5, aes(color = sig)) +
  scale_color_manual(breaks = c("cont.", "ns", "SBRT"),
                     values = c("#1e88e5", "grey", "#d81b60")) +
  guides(color = guide_legend(override.aes = list(size = 3,
                                                  shape = 15),
                              nrow = 1)) +
  labs(color = "DGE") +
  xlab(expression("control    " %<-% "    log2FC    " %->% "    SBRT")) +
  geom_text_repel(data = filter(deg, ID %in% c("HLA-A", "CA9")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -2,
                  nudge_y = 50) +
  geom_text_repel(data = filter(deg, ID %in% c("CALR")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = -2,
                  nudge_y = 100) +
  geom_text_repel(data = filter(deg, ID %in% c("IFI27", "CFH")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 3,
                  nudge_y = -10) +
  geom_text_repel(data = filter(deg, ID %in% c("IFI16", "IFI6", "BAX")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 4,
                  nudge_y = 40) +
  geom_text_repel(data = filter(deg, ID %in% c("MUC1")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 6,
                  nudge_y = 14) +
  geom_text_repel(data = filter(deg, ID %in% c("HLA-E")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 6,
                  nudge_y = -3) +
  geom_text_repel(data = filter(deg, ID %in% c("BAD")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 4,
                  nudge_y = -5) +
  geom_text_repel(data = filter(deg, ID %in% c("IFNGR1")),
                  aes(label = ID), fontface = "italic",
                  min.segment.length = 0,
                  size = 2.25,
                  nudge_x = 4,
                  nudge_y = 5) +
  scale_x_continuous(limits = c(-6, 9),
                     breaks = seq(-9, 9, 3)) +
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
  filter(grepl("gen proc|Cross-pres|DNA D|DNA Repair|gamma sig|Apop", Description))
sub_GSEA$Description <- factor(sub_GSEA$Description,
                               levels = sub_GSEA$Description[order(sub_GSEA$NES, decreasing = F)])

########## tumor GSEA ##########
png(filename = "Fig7_gsea.png", width = 850, height = 750, res = 300)

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

### probe lists
res <- c("MX1", "IFI44", "CA2", "MCL1", "TIMP3", "LAMP3", "LGALS3BP", "IFI27", "OAS1", "STAT1", "IFIT1", "TRIM14", "IRF7", "USP18", "BST2", "CXCL10", "DDX60", "HERC6", "HLA-B", "HLA-C", "IFI35", "IFI44L", "ISG15", "LY6E", "MX2", "OAS3", "OASL", "PLSCR1", "HSD17B1", "CCNA1", "CXCL1", "GALC", "IFI6", "ROBO1", "SLC6A15", "THBS1")
res <- AnnotationDbi::select(hs, keys = res, columns = c("ENTREZID"), keytype = "SYMBOL") # match gene symbol to ENTREZID
res <- res[which(is.na(res$ENTREZID) == F), "ENTREZID"] # filter genes with missing ENTREZID

hal <- c("TNFRSF10", "IRF9", "PSMB9", "EPSTI1", "PARP12", "TRIM25", "LAP3", "CASP7", "UPP1", "B2M", "IRF4", "SRI", "NFKB1A", "IFIT2", "OAS2", "TAP1","EIF2AK2", "RSAD2", "IRF1", "XAF1", "SP110", "PSMB8", "IFITM3", "GBP4", "IRF8", "PML", "IFIH1", "UBE2L6", "ADAR", "STAT2", "CXCL9", "IL10RA", "PLA2G4A", "TRIM21", "PTGS2", "C1S", "DDX58", "IL15", "NLRC5", "NMI", "IDO1", "PSMB10", "CXCL11", "ITGB7", "SAMHD1", "CMPK2", "SAMD9L", "RTP4", "PTPN2", "PARP14", "TNFAIP2", "IFITM2", "SOCS1", "CASP1", "ICAM1", "WARS", "PSME1", "ISG20", "IRF2", "FCGR1A", "MARCH1", "SOCS3", "JAK2", "HLA-DMA", "TNFAIP6", "TRIM26", "VCAM1", "CD274", "CIITA", "CIITA", "NAMPT", "SELP", "GPR18", "FPR1", "PRIC285", "PSME2", "SERPING1", "CCL5", "RNF31", "SOD2", "PSMA3", "RNF213", "PELI1", "CFB", "CD86", "TXNIP", "HLA-DQA1", "GCH1", "PNP", "CCL7", "PTPN6", "SPPL2A", "IL4R","PNPT1", "DHX58" ,"BTG1", "CASP8", "IFI30", "CCL2", "FGL2", "SECTM1", "IL15RA", "CD40", "TRAFD1", "HLA-DRB1" ,"GBP6", "LCP2", "MT2A", "RIPK1", "KLRK1", "PSMB2", "TDRD7", "HIF1A", "EIF4E3", "VAMP8", "PFKP", "CD38", "ZBP1", "BANK1", "TOR1B", "RBCK1","PDE4B", "MVP", "IL7", "BPGM", "FTSJD2", "AUTS2", "RIPK", "CD69", "MYD88", "PSMA2", "PIM1", "NOD1", "CFH", "TAPBP", "SLC25A28", "PTPN1", "TNFAIP3", "SSPN", "NUP93","MTHFD2", "CDKN1A", "NFKB1", "BATF2", "LATS2", "IRF5", "SLAMF7", "ISOC1", "P2RY14", "STAT3", "NCOA3", "HLA-A", "IL6", "GZMA", "IFNAR2", "CD74", "RAPGEF6", "GASP4", "FAS", "OGFR", "ARL4A", "LYSMD2", "CSF2RB", "ST3GAL5", "C1R", "CASP3", "CMKLR1", "METTL7B", "ST8SIA4", "XCL1", "IL2RB", "VAMP5", "IL18BP", "ZNFX1", "ARID5B", "APOL6", "STAT4")
hal <- AnnotationDbi::select(hs, keys = hal, columns = c("ENTREZID"), keytype = "SYMBOL") # match gene symbol to ENTREZID
hal <- hal[which(is.na(hal$ENTREZID) == F), "ENTREZID"] # filter genes with missing ENTREZID

### calculate GSEA of lists
score_lists <- list("IFNG hallmark" = hal,
                    "IFNG resistance" = res) # rename "Pathway" to name of choice

fgsea(pathways = score_lists,
      stats    = ent_ord,
      minSize  = 1,
      maxSize  = 500)

statsAdj <- ent_ord / max(abs(ent_ord))

########## IFNG hallmark GSEA ##########
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

png(filename = "SFig8_ih.png", width = 750, height = 350, res = 300)

ggplot(toPlot, aes(x = x, y = y)) +
  xlab(expression("SBRT    " %<-% "    Gene rank    " %->% "    control")) +
  ylab("Enrichmen\nscore") +
  ggtitle("Hallmark IFNG, NES = 1.26, padj < 0.0005") +
  geom_point(color = "green", size = 0.1) +
  geom_hline(yintercept = max(tops), color = "red", linetype = "dashed") +
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") + theme_classic() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff/3,
                             xend = x, yend = diff/3),
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

########## IFNG resistance GSEA ##########
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

png(filename = "SFig8_ir.png", width = 750, height = 350, res = 300)

ggplot(toPlot, aes(x = x, y = y)) +
  xlab(expression("SBRT    " %<-% "    Gene rank    " %->% "    control")) +
  ylab("Enrichment\nscore") +
  ggtitle("IFNG resistance, NES = 1.37, padj < 0.005") +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
  geom_point(color = "green", size = 0.1) +
  geom_hline(yintercept = max(tops), color = "red", linetype = "dashed") +
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(color = "green") + theme_classic() +
  geom_segment(data = data.frame(x = pathway),
               mapping = aes(x = x, y = -diff/3,
                             xend = x, yend = diff/3),
               size = 0.25) +
  theme(axis.text.x = element_text(size = 7,
                                   color = "white"),
        axis.text.y = element_text(size = 7,
                                   color = "black"),
        axis.title.x = element_text(size = 7,
                                    color = "black"),
        axis.title.y = element_text(size = 7,
                                    color = "black"),
        plot.title = element_text(size = 7))

dev.off()

########## IFNG response score TSNE ##########
pt <- data.frame(
  cbind(x = reducedDim(tum, "TSNE")[, 1],
        y = reducedDim(tum, "TSNE")[, 2],
        colData(tum)))
hal <- c("TNFRSF10", "IRF9", "PSMB9", "EPSTI1", "PARP12", "TRIM25", "LAP3", "CASP7", "UPP1", "B2M", "IRF4", "SRI", "NFKB1A", "IFIT2", "OAS2", "TAP1","EIF2AK2", "RSAD2", "IRF1", "XAF1", "SP110", "PSMB8", "IFITM3", "GBP4", "IRF8", "PML", "IFIH1", "UBE2L6", "ADAR", "STAT2", "CXCL9", "IL10RA", "PLA2G4A", "TRIM21", "PTGS2", "C1S", "DDX58", "IL15", "NLRC5", "NMI", "IDO1", "PSMB10", "CXCL11", "ITGB7", "SAMHD1", "CMPK2", "SAMD9L", "RTP4", "PTPN2", "PARP14", "TNFAIP2", "IFITM2", "SOCS1", "CASP1", "ICAM1", "WARS", "PSME1", "ISG20", "IRF2", "FCGR1A", "MARCH1", "SOCS3", "JAK2", "HLA-DMA", "TNFAIP6", "TRIM26", "VCAM1", "CD274", "CIITA", "CIITA", "NAMPT", "SELP", "GPR18", "FPR1", "PRIC285", "PSME2", "SERPING1", "CCL5", "RNF31", "SOD2", "PSMA3", "RNF213", "PELI1", "CFB", "CD86", "TXNIP", "HLA-DQA1", "GCH1", "PNP", "CCL7", "PTPN6", "SPPL2A", "IL4R","PNPT1", "DHX58" ,"BTG1", "CASP8", "IFI30", "CCL2", "FGL2", "SECTM1", "IL15RA", "CD40", "TRAFD1", "HLA-DRB1" ,"GBP6", "LCP2", "MT2A", "RIPK1", "KLRK1", "PSMB2", "TDRD7", "HIF1A", "EIF4E3", "VAMP8", "PFKP", "CD38", "ZBP1", "BANK1", "TOR1B", "RBCK1","PDE4B", "MVP", "IL7", "BPGM", "FTSJD2", "AUTS2", "RIPK", "CD69", "MYD88", "PSMA2", "PIM1", "NOD1", "CFH", "TAPBP", "SLC25A28", "PTPN1", "TNFAIP3", "SSPN", "NUP93","MTHFD2", "CDKN1A", "NFKB1", "BATF2", "LATS2", "IRF5", "SLAMF7", "ISOC1", "P2RY14", "STAT3", "NCOA3", "HLA-A", "IL6", "GZMA", "IFNAR2", "CD74", "RAPGEF6", "GASP4", "FAS", "OGFR", "ARL4A", "LYSMD2", "CSF2RB", "ST3GAL5", "C1R", "CASP3", "CMKLR1", "METTL7B", "ST8SIA4", "XCL1", "IL2RB", "VAMP5", "IL18BP", "ZNFX1", "ARID5B", "APOL6", "STAT4")
pt$hal <- rowSums(t(logcounts(tum[rownames(tum) %in% hal, ])))
pt$hal[pt$hal > 100] <- 100
res <- c("MX1", "IFI44", "CA2", "MCL1", "TIMP3", "LAMP3", "LGALS3BP", "IFI27", "OAS1", "STAT1", "IFIT1", "TRIM14", "IRF7", "USP18", "BST2", "CXCL10", "DDX60", "HERC6", "HLA-B", "HLA-C", "IFI35", "IFI44L", "ISG15", "LY6E", "MX2", "OAS3", "OASL", "PLSCR1", "HSD17B1", "CCNA1", "CXCL1", "GALC", "IFI6", "ROBO1", "SLC6A15", "THBS1")
pt$res <- rowSums(t(logcounts(tum[rownames(tum) %in% res, ])))
pt$res[pt$res > 30] <- 30

png(filename = "SFig8_sc.png", width = 600, height = 900, res = 300)

grid.arrange(grobs = list(
ggplot(pt, aes(x, y, color = hal)) +
  geom_point(size = 0.1) +
  xlab(expression("SBRT    " %<-% "    Gene rank    " %->% "    control")) +
  ylim(c(-26, -8)) +
  xlim(c(-33, -15)) +
  scale_color_gradient2(low = "grey90", mid = "#003831", high = "#efb21e",
                        breaks = c(25, 50, 75, 100), labels = c("25", "50", "75", "> 100")) +
  guides(color = guide_colorbar(barwidth = 0.25,
                                barheight = 3,
                                title.position = "top",
                                title = "Hallmark IFNG\nscore"),
         alpha = F) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.title = element_text(size = 7, hjust = 0),
        legend.text = element_text(size = 7),
        legend.position = "right")
,
ggplot(pt, aes(x, y, color = res)) +
  geom_point(size = 0.1) +
  ylim(c(-26, -8)) +
  xlim(c(-33, -15)) +
  scale_color_gradient2(low = "grey90", mid = "#003831", high = "#efb21e",
                        breaks = c(10, 20, 30), labels = c("10", "20", "> 30")) +
  guides(color = guide_colorbar(barwidth = 0.25,
                                barheight = 3,
                                title.position = "top",
                                title = "IFNG resistance\nscore"),
         alpha = F) +
  theme_void() +
  theme(aspect.ratio = 1,
        legend.title = element_text(size = 7, hjust = 0),
        legend.text = element_text(size = 7),
        legend.position = "right")
), nrow = 2)

dev.off()

pt$hal <- rowSums(t(logcounts(tum[rownames(tum) %in% hal, ])))
pt$res <- rowSums(t(logcounts(tum[rownames(tum) %in% res, ])))

cor.test(pt$hal, pt$res, method=c("pearson"))

png(filename = "SFig8_res_hal.png", width = 1200, height = 900, res = 300)
ggplot(pt) +
  geom_point(aes(x = hal, y = res, color = treatment), size = 0.25) +
  geom_smooth(method = "lm", aes(x = hal, y = res), color = "black") +
  labs(x = "Hallmark IFNG", y = "IFNG resistance") +
  scale_x_continuous(expand= c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("#1e88e5", "#d81b60")) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7, hjust = 0),
        legend.text = element_text(size = 7),
        legend.position = "right")
dev.off()
