library(nichenetr)
library(tidyverse)
library(KEGGREST)
library(SingleCellExperiment)
library(scater)
library(limma)
library(ComplexHeatmap)

# load data
## single cell data
#tnk <- readRDS("~/_Projects/Chow Nat Comm/data/tnk_FZ02.RDS")
#tum <- readRDS("~/_Projects/Chow Nat Comm/data/tum_FZ02.RDS")
sce <- SingleCellExperiment()
for(file_ in list.files("~/_Projects/Chow Nat Comm/data/sce_final", full.names = T)){
  if(nrow(sce) == 0){
    sce <- readRDS(file_)
  } else {
    sce <- cbind(sce, readRDS(file_))
  }
}

tnk <- sce[, sce$class == "tnk"]
tum <- sce[, sce$class == "tum"]

## ligand target gene data
ltm <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# saveRDS(ltm, "ltm.rds")
# ltm <- readRDS("ltm.rds")
## curated ligand-receptor network
lrn <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
## curated weighted network
wtn <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
# saveRDS(wtn, "wtn.rds")
# wtn <- readRDS("wtn.rds")
## KEGG gene sets of interest
# CD8 source RCC target
apop <- sub(";.*", "", keggGet("hsa04210")[[1]]$GENE[seq(2,length(keggGet("hsa04210")[[1]]$GENE), 2)]) # hsa04210 Apoptosis
agpr <- sub(";.*", "", keggGet("hsa04612")[[1]]$GENE[seq(2,length(keggGet("hsa04612")[[1]]$GENE), 2)]) # hsa04612 Antigen processing and presentation
# RCC source CD8 target
# pdcp <- sub(";.*", "", keggGet("hsa05235")[[1]]$GENE[seq(2,length(keggGet("hsa05235")[[1]]$GENE), 2)]) # hsa05235 PD-L1 expression and PD-1 checkpoint pathway in cancer
# tcrs <- sub(";.*", "", keggGet("hsa04660 ")[[1]]$GENE[seq(2,length(keggGet("hsa04660 ")[[1]]$GENE), 2)]) # hsa04660 T cell receptor signaling pathway
# filter CD8
cd8 <- tnk[, grepl("cd8", tnk$subclass)] # take only CD8 clusters
# clean up
rm(tnk)
## calculate DEG of CD8 and RCC by treatment group
# get dge w/in tumors
counts <- log2(1 + calculateCPM(tum))
design <- model.matrix(~I(treatment == "SBRT"), colData(tum))
l <- lmFit(counts, design)
e <- eBayes(l)
deg_tum <- topTable(e, n = Inf)
deg_tum$sig <- ifelse(abs(deg_tum$logFC) < 1 |
                        deg_tum$adj.P.Val > 0.05,
                      "ns",
                      ifelse(deg_tum$logFC > 1,
                             "SBRT",
                             "cont."))
rm(counts, design, l, e)

# get dge w/in CD8
counts <- log2(1 + calculateCPM(cd8))
design <- model.matrix(~I(treatment == "SBRT"), colData(cd8))
l <- lmFit(counts, design)
e <- eBayes(l)
deg_cd8 <- topTable(e, n = Inf)
deg_cd8$sig <- ifelse(abs(deg_cd8$logFC) < 1 |
                        deg_cd8$adj.P.Val > 0.05,
                      "ns",
                      ifelse(deg_cd8$logFC > 1,
                             "SBRT",
                             "cont."))
rm(counts, design, l, e)

# clean up
colData(cd8) <- colData(cd8)[, c("treatment", "sample", "subclass")]
colData(tum) <- colData(tum)[, c("treatment", "sample")]

########## CD8 source RCC target ##########
# genes of interest
tum_bg <- rownames(tum[rowMeans(counts(tum) > 0) > 0.1, ]) %>%
  .[. %in% rownames(ltm)] # background genes expressed in > 10% tumor cells and present in ltm
cd8_lg <- rownames(cd8[rowMeans(counts(cd8) > 0) > 0.1, ]) %>%
  .[. %in% c(lrn$from, colnames(ltm))] # ligand genes expressed in > 10% CD8 and present in ltm and verified ligand

# predict relevant ligands
la <- predict_ligand_activities(geneset = unique(apop, agpr), # genes of interest
                                background_expressed_genes = tum_bg, # all available target genes
                                ligand_target_matrix = ltm, # pre-calculated matrix
                                potential_ligands = cd8_lg) # potential ligands

cd8_lg_major <- la[order(la$pearson, decreasing = T), ][c(1:10, 39, 67), ] # take top predicted ligands for given targets

pearson <- cd8_lg_major

rm(la, tum_bg, cd8_lg, cd8_lg_major)

# predicted physical interactions between CD8 ligand and tum receptors
tum_rc_major <- unique(lrn[lrn$from %in% pearson$test_ligand &
                             lrn$database == "kegg", ]$to) # take kegg curated receptors for major ligands

cd8_tum_wtn <- wtn$lr_sig %>%
  filter(from %in% pearson$test_ligand &
           to %in% tum_rc_major)
cd8_tum_wtn <- cd8_tum_wtn %>% spread(to, weight)
lr_int <- as.matrix(cd8_tum_wtn[, 2:16])
rownames(lr_int) <- c(cd8_tum_wtn$from)
lr_int <- rbind(lr_int,
                GLG1 = c(rep(NA, 15)),
                TIGIT = c(rep(NA, 15)),
                GZMB = c(rep(NA, 15)))

rm(cd8_tum_wtn)

# predicted genetic interactions
top_lg_tg <- c() # get top predicted ligand gene targets
for(x in pearson$test_ligand){
  top_lg_tg <- c(top_lg_tg,
                 rownames(ltm)[order(ltm[, x], decreasing = T)[1:500]]) # top 500 predicted target
}
top_lg_tg <- unique(top_lg_tg)

genes <- c(apop, agpr)[c(apop, agpr) %in% unique(top_lg_tg)] # genes of interest are putative ligand targets
genes <- genes[genes %in% deg_tum$ID[deg_tum$sig != "ns"]] # genes of interest are deg between treatment groups
genes <- unique(genes)

lt_int <- t(ltm[genes, as.character(pearson$test_ligand)]) # subset ligands and their targets of interest

rm(top_lg_tg, x)

# ligand expression data
lg <- logcounts(cd8)[pearson$test_ligand, ]
lgdat <- data.frame(cont.avg = rowMeans(lg[, cd8$treatment == "control"]), cont.pos = rowMeans(lg[, cd8$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(lg[, cd8$treatment == "SBRT"]), sbrt.pos = rowMeans(lg[, cd8$treatment == "SBRT"] > 0))
lgdeg <- deg_cd8[deg_cd8$ID %in% rownames(lgdat), c("logFC", "adj.P.Val", "ID")]
rownames(lgdeg) <- lgdeg$ID
lgdat <- cbind(lgdat, lgdeg[rownames(lgdat), ])

rm(lg, lgdeg)

# receptor expression data
rc <- logcounts(tum)[tum_rc_major, ]
rcdat <- data.frame(cont.avg = rowMeans(rc[, tum$treatment == "control"]), cont.pos = rowMeans(rc[, tum$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(rc[, tum$treatment == "SBRT"]), sbrt.pos = rowMeans(rc[, tum$treatment == "SBRT"] > 0))
rcdeg <- deg_tum[deg_tum$ID %in% rownames(rcdat), c("logFC", "adj.P.Val", "ID")]
rownames(rcdeg) <- rcdeg$ID
rcdat <- cbind(rcdat, rcdeg[rownames(rcdat), ])

rm(rc, rcdeg)

# target gene expression data
tg <- logcounts(tum)[genes, ]
tgdat <- data.frame(cont.avg = rowMeans(tg[, tum$treatment == "control"]), cont.pos = rowMeans(tg[, tum$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(tg[, tum$treatment == "SBRT"]), sbrt.pos = rowMeans(tg[, tum$treatment == "SBRT"] > 0))
tgdeg <- deg_tum[deg_tum$ID %in% rownames(tgdat), c("logFC", "adj.P.Val", "ID")]
rownames(tgdeg) <- tgdeg$ID
tgdat <- cbind(tgdat, tgdeg[rownames(tgdat), ])

rm(tg, tgdeg)

# heatmap
lr_int <- lr_int[pearson$test_ligand,
                 colnames(lr_int)[order(colSums(lr_int, na.rm = T), decreasing = T)]]

# annotation data
lgdat$logFC[lgdat$logFC > 1] <- 1 # scale logFC
lgdat$logFC[lgdat$logFC < -1] <- -1
lgdat$sc.adj.p <- lgdat$adj.P.Val
lgdat$sc.adj.p[lgdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
lgdat$sc.adj.p <- -log10(lgdat$sc.adj.p)/max(-log10(lgdat$sc.adj.p))

lgsc <- c()
for(x in sort(unique(cd8$subclass))){
  a <- logcounts(cd8)[pearson$test_ligand, cd8$subclass == x]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  lgsc <- cbind(lgsc, a)
}
colnames(lgsc) <- c("CD8_TCF7", "CD8_ITM2C", "CD8_CCR5", "CD8_HAVCR2", "CD8_MKI67")
lgsc <- t(scale(t(lgsc)))
lgsc[lgsc > 2] <- 2
lgsc[lgsc < -2] <- -2
rm(a, x)

rcdat <- rcdat[colnames(lr_int), ]
rcdat$logFC[rcdat$logFC > 1] <- 1 # scale logFC
rcdat$logFC[rcdat$logFC < -1] <- -1
rcdat$sc.adj.p <- rcdat$adj.P.Val
rcdat$sc.adj.p[rcdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
rcdat$sc.adj.p <- -log10(rcdat$sc.adj.p)/max(-log10(rcdat$sc.adj.p))

tgdat <- tgdat[colnames(lt_int), ]
tgdat$logFC[tgdat$logFC > 1] <- 1 # scale logFC
tgdat$logFC[tgdat$logFC < -1] <- -1
tgdat$sc.adj.p <- tgdat$adj.P.Val
tgdat$sc.adj.p[tgdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
tgdat$sc.adj.p <- -log10(tgdat$sc.adj.p)/max(-log10(tgdat$sc.adj.p))

# color scales
col_lr_int <- circlize::colorRamp2(seq(0, 1, 0.2),
                                   colorRampPalette(c("darkorange", "darkorange4"))(6))
col_lt_int <- circlize::colorRamp2(seq(0, 0.008, 0.001),
                                   colorRampPalette(c("white", "green4"))(9))
col_deg <- circlize::colorRamp2(seq(-1, 1, 0.5),
                                colorRampPalette(c("#1e88e5", "white", "#d81b60"))(5))
col_exp <- circlize::colorRamp2(seq(0, 5, 0.5),
                                colorRampPalette(c("grey90", "#003831"))(11))
col_fun <- circlize::colorRamp2(c(-2, -1, 0, 1, 2),
                                colorRampPalette(c("purple", "black", "yellow"))(5))

lg_ha <- rowAnnotation(
  cont_expression = anno_simple(rep(0, 12),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(lgdat$cont.avg)),
                                pt_size = unit(lgdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 12),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(lgdat$sbrt.avg)),
                                pt_size = unit(lgdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 12),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(lgdat$logFC)),
                    pt_size = unit(lgdat$sc.adj.p * 3, "mm")),
  subclass = lgsc,
  PCC = pearson$pearson,
  col = list(
    PCC = circlize::colorRamp2(seq(0.1, 0.35, 0.05),
                               colorRampPalette(c("white", "midnightblue"))(6)),
    subclass = col_fun),
  gp = gpar(col = "black"),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 7),
  show_legend = F
)

rc_ha <- HeatmapAnnotation(
  cont_expression = anno_simple(rep(0, 15),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(rcdat$cont.avg)),
                                pt_size = unit(rcdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 15),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(rcdat$sbrt.avg)),
                                pt_size = unit(rcdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 15),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(rcdat$logFC)),
                    pt_size = unit(rcdat$sc.adj.p * 3, "mm")),
  annotation_label = c("", "", ""),
  annotation_name_gp = gpar(fontsize = 7)
)

tg_ha <- HeatmapAnnotation(
  cont_expression = anno_simple(rep(0, 30),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(tgdat$cont.avg)),
                                pt_size = unit(tgdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 30),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(tgdat$sbrt.avg)),
                                pt_size = unit(tgdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 30),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(tgdat$logFC)),
                    pt_size = unit(tgdat$sc.adj.p * 3, "mm")),
  annotation_name_gp = gpar(fontsize = 7)
)

set.seed(415); cd8_rcc_lr <- Heatmap(lr_int,
                       na_col = "grey80",
                       show_heatmap_legend = F,
                       col = col_lr_int,
                       left_annotation = lg_ha,
                       top_annotation = rc_ha,
                       cluster_rows = F,
                       cluster_columns = F,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title = "CD8 ligands",
                       row_title_gp = gpar(fontsize = 7),
                       column_title = "RCC receptors",
                       column_title_gp = gpar(fontsize = 7),
                       width = ncol(lr_int)*unit(3, "mm"), 
                       height = nrow(lr_int)*unit(3, "mm"),
                       rect_gp = gpar(col = "black"))

set.seed(415); cd8_rcc_lt <- Heatmap(lt_int,
                       na_col = "white",
                       show_heatmap_legend = F,
                       col = col_lt_int,
                       top_annotation = tg_ha,
                       cluster_rows = F,
                       cluster_columns = F,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       row_title = "CD8 ligands",
                       row_title_gp = gpar(fontsize = 7),
                       column_title = "RCC target genes",
                       column_title_gp = gpar(fontsize = 7),
                       width = ncol(lt_int)*unit(3, "mm"), 
                       height = nrow(lt_int)*unit(3, "mm"),
                       rect_gp = gpar(col = "black"))

png(filename = "Fig7_ip.png", width = 2550, height = 1000, res = 300)
cd8_rcc_lr + cd8_rcc_lt
dev.off()

p_lgd = Legend(col_fun = circlize::colorRamp2(seq(0.1, 0.35, 0.05),
                                              colorRampPalette(c("white", "midnightblue"))(6)),
               title = "Pearson corellation coefficient",
               at = seq(0.1, 0.35, 0.05),
               direction = "horizontal",
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               grid_height = unit(0.25, "cm"),
               legend_width = unit(3, "cm"))
lr_lgd = Legend(col_fun = col_lr_int,
               title = "Interaction potential",
               at = seq(0, 1, 0.2),
               direction = "horizontal",
               title_gp = gpar(fontsize = 7),
               labels_gp = gpar(fontsize = 7),
               grid_height = unit(0.25, "cm"),
               legend_width = unit(3, "cm"))
lt_lgd = Legend(col_fun = col_lt_int,
                title = "Regulatory potential",
                at = seq(0, 0.008, 0.002),
                direction = "horizontal",
                title_gp = gpar(fontsize = 7),
                labels_gp = gpar(fontsize = 7),
                grid_height = unit(0.25, "cm"),
                legend_width = unit(3, "cm"))
sc_lgd = Legend(col_fun = col_fun,
                title = "Z score",
                at = seq(-2, 2, 1),
                labels = c("< -2", "-1", "0", "1", "> 2"),
                direction = "horizontal",
                title_gp = gpar(fontsize = 7),
                labels_gp = gpar(fontsize = 7),
                grid_height = unit(0.25, "cm"),
                legend_width = unit(3, "cm"))

png(filename = "Fig7_ip_lg1.png", width = 500, height = 600, res = 300)
draw(packLegend(sc_lgd, p_lgd, lr_lgd, lt_lgd))
dev.off()

exp_lgd = Legend(col_fun = col_exp,
                 title = "Mean\nlognormcounts",
                 at = seq(0, 5, 1),
                 direction = "vertical",
                 title_gp = gpar(fontsize = 7),
                 labels_gp = gpar(fontsize = 7),
                 grid_width = unit(0.25, "cm"),
                 grid_height = unit(0.25, "cm"))
exp_pos_lgd = Legend(title = "\n% positive",
                     labels = paste0(c(0, 33, 66, 100), "%"),
                     type = "points",
                     title_gp = gpar(fontsize = 7),
                     labels_gp = gpar(fontsize = 7),
                     pch = 21,
                     size = unit(0:3, "mm"))
deg_lgd = Legend(col_fun = col_deg,
                 title = "\nDGE, log2FC",
                 at = seq(-1, 1, 0.5),
                 labels = c("< -1", "-0.5", "0", "0.5", "> 1"),
                 direction = "vertical",
                 title_gp = gpar(fontsize = 7),
                 labels_gp = gpar(fontsize = 7),
                 grid_width = unit(0.25, "cm"),
                 grid_height = unit(0.33, "cm"))
deg_pos_lgd = Legend(title = "\nSignificance\n-log(adj.p.value)",
                     labels = c("p = 0.25", "p < 0.05"),
                     type = "points",
                     title_gp = gpar(fontsize = 7),
                     labels_gp = gpar(fontsize = 7),
                     pch = 21,
                     size = unit(c(1.5, 3), "mm"))

png(filename = "Fig7_ip_lg2.png", width = 500, height = 600, res = 300)
draw(packLegend(exp_lgd, exp_pos_lgd, deg_lgd, deg_pos_lgd,
                direction = "horizontal",
                column_gap = unit(5, "mm"),
                max_width = unit(5, "cm")))
dev.off()

########## CD8 source Mp target ##########
mye <- readRDS("~/_Projects/Chow Nat Comm/data/mye_FZ02.RDS")
mac <- mye[, grepl("mac", mye$subclass)]
rm(mye)

# get dge w/in Mp
counts <- log2(1 + calculateCPM(mac))
design <- model.matrix(~I(treatment == "SBRT"), colData(mac))
l <- lmFit(counts, design)
e <- eBayes(l)
deg_mac <- topTable(e, n = Inf)
deg_mac$sig <- ifelse(abs(deg_mac$logFC) < 1 |
                        deg_mac$adj.P.Val > 0.05,
                      "ns",
                      ifelse(deg_mac$logFC > 1,
                             "SBRT",
                             "cont."))

# genes of interest
mac_bg <- rownames(mac[rowMeans(counts(mac) > 0) > 0.1, ]) %>%
  .[. %in% rownames(ltm)] # background genes expressed in > 10% macrophages and present in ltm
cd8_lg <- rownames(cd8[rowMeans(counts(cd8) > 0) > 0.1, ]) %>%
  .[. %in% colnames(ltm)] # ligand genes expressed in > 10% CD8 and present in ltm

# predict relevant ligands
mactar <- c("JAK1", "JAK2", "STAT1",
            "NOS1", "NOS3", "IL12A", # KEGG 05140
            "IRF8", "IRF1", # Dror Molecular Immunology 2007
            "SOCS3", "SOCS7", "MAP4K4", 
            "NFKB1", "NFKB2", "IL1B", "TNF", "IL6",
            "MMP9", "LIF", "CXCL16")

la <- predict_ligand_activities(geneset = mactar,
                                background_expressed_genes = mac_bg, # all available target genes
                                ligand_target_matrix = ltm, # pre-calculated matrix
                                potential_ligands = c(cd8_lg, "IL10", "IL15")) # potential ligands

cd8_lg_major <- la[order(la$pearson, decreasing = T), ][c(1:5, 7, 20, 31), ] # take top ligands and cytokines from Braun Cancer Cell 2021

pearson <- cd8_lg_major

rm(la, mac_bg, cd8_lg, cd8_lg_major)

# predicted physical interactions between CD8 ligand and Mp receptors
mac_rc_major <- unique(lrn[lrn$from %in% pearson$test_ligand &
                             lrn$database == "kegg", ]$to) # take kegg curated receptors for major ligands

cd8_mac_wtn <- wtn$lr_sig %>%
  filter(from %in% pearson$test_ligand &
           to %in% mac_rc_major)
cd8_mac_wtn <- cd8_mac_wtn %>% spread(to, weight)
lr_int <- as.matrix(cd8_mac_wtn[, 2:13])
rownames(lr_int) <- c(cd8_mac_wtn$from)
lr_int <- lr_int[pearson$test_ligand, c("IFNGR1", "IFNGR2",
                                        "TNFRSF1A", "TNFRSF1B",
                                        "CCR1", "CCR3", "CCR5",
                                        "IL15RA",
                                        "IL2RB", "IL2RG",
                                        "IL10RA", "IL10RB")]

# predicted genetic interactions
top_lg_tg <- c() # get top predicted ligand gene targets
for(x in pearson$test_ligand){
  top_lg_tg <- c(top_lg_tg,
                 rownames(ltm)[order(ltm[, x], decreasing = T)[1:500]]) # top 500 predicted target
}
top_lg_tg <- unique(top_lg_tg)

lt_int <- t(ltm[mactar[mactar %in% rownames(ltm)], as.character(pearson$test_ligand)]) # subset ligands and their targets of interest

rm(top_lg_tg, x)

# ligand expression data
lg <- logcounts(cd8)[pearson$test_ligand, ]
lgdat <- data.frame(cont.avg = rowMeans(lg[, cd8$treatment == "control"]), cont.pos = rowMeans(lg[, cd8$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(lg[, cd8$treatment == "SBRT"]), sbrt.pos = rowMeans(lg[, cd8$treatment == "SBRT"] > 0))
lgdeg <- deg_cd8[deg_cd8$ID %in% rownames(lgdat), c("logFC", "adj.P.Val", "ID")]
rownames(lgdeg) <- lgdeg$ID
lgdat <- cbind(lgdat, lgdeg[rownames(lgdat), ])

rm(lg, lgdeg)

# receptor expression data
rc <- logcounts(mac)[mac_rc_major, ]
rcdat <- data.frame(cont.avg = rowMeans(rc[, mac$treatment == "control"]), cont.pos = rowMeans(rc[, mac$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(rc[, mac$treatment == "SBRT"]), sbrt.pos = rowMeans(rc[, mac$treatment == "SBRT"] > 0))
rcdeg <- deg_mac[deg_mac$ID %in% rownames(rcdat), c("logFC", "adj.P.Val", "ID")]
rownames(rcdeg) <- rcdeg$ID
rcdat <- cbind(rcdat, rcdeg[rownames(rcdat), ])

rm(rc, rcdeg)

# target gene expression data
tg <- logcounts(mac)[colnames(lt_int), ]
tgdat <- data.frame(cont.avg = rowMeans(tg[, mac$treatment == "control"]), cont.pos = rowMeans(tg[, mac$treatment == "control"] > 0),
                    sbrt.avg = rowMeans(tg[, mac$treatment == "SBRT"]), sbrt.pos = rowMeans(tg[, mac$treatment == "SBRT"] > 0))
tgdeg <- deg_mac[deg_mac$ID %in% rownames(tgdat), c("logFC", "adj.P.Val", "ID")]
rownames(tgdeg) <- tgdeg$ID
tgdat <- cbind(tgdat, tgdeg[rownames(tgdat), ])

rm(tg, tgdeg)

# heatmap
# annotation data
lgdat$logFC[lgdat$logFC > 1] <- 1 # scale logFC
lgdat$logFC[lgdat$logFC < -1] <- -1
lgdat$sc.adj.p <- lgdat$adj.P.Val
lgdat$sc.adj.p[lgdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
lgdat$sc.adj.p <- -log10(lgdat$sc.adj.p)/max(-log10(lgdat$sc.adj.p))

rcdat <- rcdat[colnames(lr_int), ]
rcdat$logFC[rcdat$logFC > 1] <- 1 # scale logFC
rcdat$logFC[rcdat$logFC < -1] <- -1
rcdat$sc.adj.p <- rcdat$adj.P.Val
rcdat$sc.adj.p[rcdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
rcdat$sc.adj.p <- -log10(rcdat$sc.adj.p)/max(-log10(rcdat$sc.adj.p))

tgdat <- tgdat[colnames(lt_int), ]
tgdat$logFC[tgdat$logFC > 1] <- 1 # scale logFC
tgdat$logFC[tgdat$logFC < -1] <- -1
tgdat$sc.adj.p <- tgdat$adj.P.Val
tgdat$sc.adj.p[tgdat$sc.adj.p < 0.05] <- 0.05 # scale adj.P
tgdat$sc.adj.p <- -log10(tgdat$sc.adj.p)/max(-log10(tgdat$sc.adj.p))

# subclass gene expression
lgsc <- c()
for(x in sort(unique(cd8$subclass))){
  a <- logcounts(cd8)[pearson$test_ligand, cd8$subclass == x]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  lgsc <- cbind(lgsc, a)
}
colnames(lgsc) <- c("CD8_TCF7", "CD8_ITM2C", "CD8_CCR5", "CD8_HAVCR2", "CD8_MKI67")
lgsc <- t(scale(t(lgsc)))
lgsc[lgsc > 2] <- 2
lgsc[lgsc < -2] <- -2

rcsc <- c()
for(x in sort(unique(mac$subclass))){
  a <- logcounts(mac)[rcdat$ID, mac$subclass == x]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  rcsc <- cbind(rcsc, a)
}
colnames(rcsc) <- c("Mp_INFG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1", "Mp_CD163", "Mp_SELENOP")
rcsc <- t(scale(t(rcsc)))
rcsc[rcsc > 2] <- 2
rcsc[rcsc < -2] <- -2

tgsc <- c()
for(x in sort(unique(mac$subclass))){
  a <- logcounts(mac)[tgdat$ID, mac$subclass == x]
  a <- cbind(rowMeans(a))
  colnames(a) <- x
  tgsc <- cbind(tgsc, a)
}
colnames(tgsc) <- c("Mp_INFG", "Mp_VEGFA", "Mp_ACP5", "Mp_CD80", "Mp_SLC40A1", "Mp_CD163", "Mp_SELENOP")
tgsc <- t(scale(t(tgsc)))
tgsc[tgsc > 2] <- 2
tgsc[tgsc < -2] <- -2

col_lr_int <- circlize::colorRamp2(seq(0, 1, 0.2),
                                   colorRampPalette(c("darkorange", "darkorange4"))(6))
col_lt_int <- circlize::colorRamp2(seq(0, 0.008, 0.001),
                                   colorRampPalette(c("white", "springgreen4"))(9))

lg_ha <- rowAnnotation(
  cont_expression = anno_simple(rep(0, 8),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(lgdat$cont.avg)),
                                pt_size = unit(lgdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 8),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(lgdat$sbrt.avg)),
                                pt_size = unit(lgdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 8),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(lgdat$logFC)),
                    pt_size = unit(lgdat$sc.adj.p * 3, "mm")),
  subclass = lgsc,
  PCC = pearson$pearson,
  col = list(
    PCC = circlize::colorRamp2(seq(0.1, 0.35, 0.05),
                               colorRampPalette(c("white", "midnightblue"))(6)),
    subclass = col_fun),
  gp = gpar(col = "black"),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 7),
  show_legend = F
)

rc_ha <- HeatmapAnnotation(
  cont_expression = anno_simple(rep(0, 12),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(rcdat$cont.avg)),
                                pt_size = unit(rcdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 12),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(rcdat$sbrt.avg)),
                                pt_size = unit(rcdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 12),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(rcdat$logFC)),
                    pt_size = unit(rcdat$sc.adj.p * 3, "mm")),
  subclass = rcsc,
  col = list(
    subclass = col_fun),
  gp = gpar(col = "black"),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 7),
  show_legend = F,
  show_annotation_name = F
)

tg_ha <- HeatmapAnnotation(
  cont_expression = anno_simple(rep(0, 19),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(tgdat$cont.avg)),
                                pt_size = unit(tgdat$cont.pos * 3, "mm")),
  SBRT_expression = anno_simple(rep(0, 19),
                                col = col_deg,
                                pch = 21,
                                pt_gp = gpar(fill = col_exp(tgdat$sbrt.avg)),
                                pt_size = unit(tgdat$sbrt.pos * 3, "mm")),
  DGE = anno_simple(rep(0, 19),
                    col = col_deg,
                    pch = 21,
                    pt_gp = gpar(fill = col_deg(tgdat$logFC)),
                    pt_size = unit(tgdat$sc.adj.p * 3, "mm")),
  subclass = tgsc,
  col = list(
    subclass = col_fun),
  gp = gpar(col = "black"),
  simple_anno_size = unit(3, "mm"),
  annotation_name_gp = gpar(fontsize = 7),
  show_legend = F
)

set.seed(415); cd8_mac_lr <- Heatmap(lr_int,
                                     na_col = "grey80",
                                     show_heatmap_legend = F,
                                     col = col_lr_int,
                                     left_annotation = lg_ha,
                                     top_annotation = rc_ha,
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     row_names_side = "left",
                                     row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                                     column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                                     row_title = "CD8 ligands",
                                     row_title_gp = gpar(fontsize = 7),
                                     column_title = "Mp receptors",
                                     column_title_gp = gpar(fontsize = 7),
                                     width = ncol(lr_int)*unit(3, "mm"), 
                                     height = nrow(lr_int)*unit(3, "mm"),
                                     rect_gp = gpar(col = "black"))

set.seed(415); cd8_mac_lt <- Heatmap(lt_int,
                                     na_col = "white",
                                     show_heatmap_legend = F,
                                     col = col_lt_int,
                                     top_annotation = tg_ha,
                                     cluster_rows = F,
                                     cluster_columns = F,
                                     row_names_side = "left",
                                     row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                                     column_names_gp = gpar(fontsize = 7, fontface = "italic"),
                                     row_title = "CD8 ligands",
                                     row_title_gp = gpar(fontsize = 7),
                                     column_title = "Mp target genes",
                                     column_title_gp = gpar(fontsize = 7),
                                     width = ncol(lt_int)*unit(3, "mm"), 
                                     height = nrow(lt_int)*unit(3, "mm"),
                                     rect_gp = gpar(col = "black"))

png(filename = "SFig7_ip.png", width = 2200, height = 1250, res = 300)
cd8_mac_lr + cd8_mac_lt
dev.off()