library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(ComplexHeatmap)

# KIRC
TCGAbiolinks:::getProjectSummary("TCGA-KIRC") # look up ccRCC TCGA

query_TCGA <- GDCquery( # evaluate available data
  project = "TCGA-KIRC", # clear cell RCC
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq", # with RNAseq data
  workflow.type = "STAR - Counts", # HTSeq estimated counts
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor") # primary tumors only
GDCdownload(query = query_TCGA) # download data of interest
kirc_tcga <- GDCprepare(query_TCGA) # unpack data
rm(query_TCGA) # clean up

saveRDS(kirc_tcga, "~/_Projects/Chow Nat Comm/tcga/kirc_tcga.RDS") # freeze 2023_0301
tcga_data <- readRDS("~/_Projects/Chow Nat Comm/data/tcga/kirc_tcga.RDS")

# multi-variate all stages
clin <- tcga_data
clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up", "ajcc_pathologic_stage")]
clin_df$stage <- clin_df$ajcc_pathologic_stage
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)
trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- ifelse(trail < 0.30, "low", "hi")

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)
ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- ifelse(ifng < 0.30, "low", "hi")

# TRAIL IFNG double
clin_df$trail_ifng_hi <- clin_df$trail.rank == "hi" & clin_df$ifng.rank == "hi"

png("SFig_strat_st_dub.png", width = 1500, height = 1000, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ stage + trail_ifng_hi, data = clin_df),
           palette = c("green4", "green1", "blue4", "blue1", "yellow4", "yellow1", "red4", "red1"),
           censor.size = 3,
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 4))
dev.off()

# stage 1-3 only
clin <- tcga_data[, tcga_data$ajcc_pathologic_stage %in% paste("Stage", c("I", "II", "III"))]

clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up", "ajcc_pathologic_stage", "gender")]
clin_df$sex <- clin_df$gender
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)

trail_ <- trail[, 1]
trail_ <- rank(trail_)
trail_ <- trail_/max(trail_)
clin_df$TNFRSF10A <- ifelse(trail_ < 0.30, "low", "hi")

trail_ <- trail[, 2]
trail_ <- rank(trail_)
trail_ <- trail_/max(trail_)
clin_df$TNFSF10 <- ifelse(trail_ < 0.30, "low", "hi")

trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- factor(ifelse(trail < 0.30, "low", "hi"), levels = c("low", "hi"))

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)

ifng_ <- ifng[, 1]
ifng_ <- rank(ifng_)
ifng_ <- ifng_/max(ifng_)
clin_df$IFNGR1 <- ifelse(ifng_ < 0.30, "low", "hi")

ifng_ <- ifng[, 2]
ifng_ <- rank(ifng_)
ifng_ <- ifng_/max(ifng_)
clin_df$IFNG <- ifelse(ifng_ < 0.30, "low", "hi")

ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- factor(ifelse(ifng < 0.30, "low", "hi"), levels = c("low", "hi"))

# TRAIL IFNG double
clin_df$trail_ifng_hi <- clin_df$trail.rank == "hi" & clin_df$ifng.rank == "hi"

png("SFig_surv_quad.png", width = 1450, height = 750, res = 300)
ggplot(data.frame(cbind(trail = clin_df$trail, ifng = clin_df$ifng)), aes(trail, ifng)) +
  geom_point(size = 0.5,
             aes(color = as.factor(paste(clin_df$ifng.rank, clin_df$trail.rank)))) +
  scale_color_manual(values = c("green3", "blue3", "darkgoldenrod3", "slategray4"),
                     labels = c("IFNG/R hi, TRAIL/R hi - 56%",
                                "IFNG/R hi, TRAIL/R lo - 14%",
                                "IFNG/R lo, TRAIL/R hi - 14%",
                                "IFNG/R lo, TRAIL/R lo - 16%")) +
  geom_vline(xintercept = 0.3, color = "red") +
  geom_hline(yintercept = 0.3, color = "red") +
  labs(color = "",
       x = "TRAIL/R",
       y = "IFNG/R") +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme_classic() +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 7,
                                 color = "black"),
        
        axis.title = element_text(size = 7,
                                 color = "black"),
        legend.text = element_text(size = 7,
                                   color = "black"))
dev.off()

# IFNG ligand
png("SFig_surv_ifngl.png", width = 500, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ IFNG, data = clin_df),
           palette = c("black", "gray60"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

coxph(Surv(overall_survival, deceased) ~ IFNG, data = clin_df) %>% 
  tbl_regression(exp = TRUE)

# IFNG receptor
png("SFig_surv_ifngr.png", width = 500, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ IFNGR1, data = clin_df),
           palette = c("black", "gray60"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

# IFNG/R score
png("SFig_surv_ifng.png", width = 1200, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ ifng.rank, data = clin_df),
           palette = c("blue1", "gray50"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T)
dev.off()
tbl_regression(coxph(Surv(overall_survival, deceased) ~ ifng.rank, data = clin_df), exp = T)

# TRAIL ligand
png("SFig_surv_traill.png", width = 500, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ TNFSF10, data = clin_df),
           palette = c("black", "gray60"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

# TRAIL receptor
png("SFig_surv_trailr.png", width = 500, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ TNFRSF10A, data = clin_df),
           palette = c("black", "gray60"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

# TRAIL/R score
png("SFig_surv_trail.png", width = 1200, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank, data = clin_df),
           palette = c("goldenrod1", "gray50"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           legend.labs = c("TRAIL/R hi", "TRAIL/R lo"),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T)
dev.off()
tbl_regression(coxph(Surv(overall_survival, deceased) ~ trail.rank, data = clin_df), exp = T)

# IFNG/R TRAIL/R score
png("SFig_surv_ifng_trail.png", width = 1200, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ ifng.rank + trail.rank,
                   data = clin_df),
           palette = c("green3", "blue3", "darkgoldenrod3", "slategray4"),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend.title = c(""),
           legend.labs = c("IFNG/R hi\nTRAIL/R hi",
                           "IFNG/R hi\nTRAIL/R lo",
                           "IFNG/R lo\nTRAIL/R hi",
                           "IFNG/R lo\nTRAIL/R lo"),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T)
dev.off()

tbl_regression(coxph(Surv(overall_survival, deceased) ~ ifng.rank + trail.rank, data = clin_df), exp = T)

# IFNG/R TRAIL/R double hi
png("Fig6_surv_dub.png", width = 1000, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df),
           palette = c("black", "green3"),
           censor.size = 3,
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend = "none",
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T)
dev.off()

tbl_regression(coxph(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df), exp = T)

lgd <- Legend(labels = c("IFNG/R TRAIL/R hi", "other"),
                title = "Strata",
                legend_gp = gpar(fill = c("green3", "black")),
                title_gp = gpar(fontsize = 7),
                title_position = "leftcenter",
                labels_gp = gpar(fontsize = 7),
                grid_height = unit(0.25, "cm"),
                grid_width = unit(0.25, "cm"),
                nrow = 2)

png(filename = "Fig6_surv_lgd.png", width = 1450, height = 250, res = 300)

draw(lgd)

dev.off()

# strata by sex
png("Fig6_surv_dub_sex.png", width = 1250, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail_ifng_hi + sex, data = clin_df),
           palette = c("goldenrod4", "green4", "goldenrod2", "green2"),
           censor.size = 3,
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  guides(colour = guide_legend(nrow = 2))
dev.off()

# LUAD stage 1-3
TCGAbiolinks:::getProjectSummary("TCGA-LUAD") # look up LUAD TCGA

query_TCGA <- GDCquery( # evaluate available data
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq", # with RNAseq data
  workflow.type = "STAR - Counts", # HTSeq estimated counts
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor") # primary tumors only
GDCdownload(query = query_TCGA) # download data of interest
luad_tcga <- GDCprepare(query_TCGA) # unpack data
rm(query_TCGA) # clean up

saveRDS(luad_tcga, "~/_Projects/Chow Nat Comm/data/tcga/luad_tcga.RDS") # freeze 2022_0617
tcga_data <- readRDS("~/_Projects/Chow Nat Comm/data/tcga/luad_tcga.RDS")

tcga_data <- tcga_data[, !is.na(tcga_data$ajcc_pathologic_stage)] # keep only data for where staging was determined
clin <- tcga_data

clin <- tcga_data[, tcga_data$ajcc_pathologic_stage != "Stage IV"] # stage 1-3 only

clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up", "ajcc_pathologic_stage")]
clin_df$stage <- str_replace(clin_df$ajcc_pathologic_stage, c("A|B"), "")
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)
trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- ifelse(trail < 0.30, "low", "hi")

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)
ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- ifelse(ifng < 0.30, "low", "hi")

png("Fig6_surv_dub_luad.png", width = 1000, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df),
           palette = c("black", "green3"),
           censor.size = 3,
           font.title = c(10),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend = "none",
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  ggtitle("LUAD")
dev.off()

# SKCM stage 1-3
TCGAbiolinks:::getProjectSummary("TCGA-SKCM") # look up SKCM TCGA

query_TCGA <- GDCquery( # evaluate available data
  project = "TCGA-SKCM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq", # with RNAseq data
  workflow.type = "STAR - Counts", # HTSeq estimated counts
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor") # primary tumors only
GDCdownload(query = query_TCGA) # download data of interest
skcm_tcga <- GDCprepare(query_TCGA) # unpack data
rm(query_TCGA) # clean up

saveRDS(skcm_tcga, "~/_Projects/Chow Nat Comm/data/tcga/skcm_tcga.RDS") # freeze 2022_0617
tcga_data <- readRDS("~/_Projects/Chow Nat Comm/data/tcga/skcm_tcga.RDS")

tcga_data <- tcga_data[, !is.na(tcga_data$ajcc_pathologic_stage)] # keep only data for where staging was determined
clin <- tcga_data[, tcga_data$ajcc_pathologic_stage != "Stage IV"] # stage 1-3 only

clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up")]
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)
trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- ifelse(trail < 0.30, "low", "hi")

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)
ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- ifelse(ifng < 0.30, "low", "hi")

png("Fig6_surv_dub_skcm.png", width = 1000, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df),
           palette = c("black", "green3"),
           censor.size = 3,
           font.title = c(10),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend = "none",
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  ggtitle("SKCM")
dev.off()

# COAD stage 1-3
TCGAbiolinks:::getProjectSummary("TCGA-COAD") # look up COAD TCGA

query_TCGA <- GDCquery( # evaluate available data
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq", # with RNAseq data
  workflow.type = "STAR - Counts", # HTSeq estimated counts
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor") # primary tumors only
GDCdownload(query = query_TCGA) # download data of interest
coad_tcga <- GDCprepare(query_TCGA) # unpack data
rm(query_TCGA) # clean up

saveRDS(skcm_tcga, "~/_Projects/Chow Nat Comm/data/tcga/coad_tcga.RDS") # freeze 2022_0628
tcga_data <- readRDS("~/_Projects/Chow Nat Comm/data/tcga/coad_tcga.RDS")

tcga_data <- tcga_data[, !is.na(tcga_data$ajcc_pathologic_stage)] # keep only data for where staging was determined
clin <- tcga_data[, !grepl("IV", tcga_data$ajcc_pathologic_stage)] # stage 1-3 only

clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up")]
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)
trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- ifelse(trail < 0.30, "low", "hi")

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)
ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- ifelse(ifng < 0.30, "low", "hi")

png("Fig6_surv_dub_coad.png", width = 1000, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df),
           palette = c("black", "green3"),
           censor.size = 3,
           font.title = c(10),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend = "none",
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  ggtitle("COAD")
dev.off()

# GBM
TCGAbiolinks:::getProjectSummary("TCGA-GBM") # look up GBM TCGA

query_TCGA <- GDCquery( # evaluate available data
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq", # with RNAseq data
  workflow.type = "STAR - Counts", # HTSeq estimated counts
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor") # primary tumors only
GDCdownload(query = query_TCGA) # download data of interest
gbm_tcga <- GDCprepare(query_TCGA) # unpack data
rm(query_TCGA) # clean up

saveRDS(skcm_tcga, "~/_Projects/Chow Nat Comm/data/tcga/gbm_tcga.RDS") # freeze 2022_0617
tcga_data <- readRDS("~/_Projects/Chow Nat Comm/data/tcga/gbm_tcga.RDS")

clin <- tcga_data # staging data unavailable

clin_df <- colData(clin)
clin_df <- clin_df[c("patient", "vital_status", "days_to_death",
                     "days_to_last_follow_up")]
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# TRAIL score, scale counts, sum counts, rank
trail <- assay(clin)[rowRanges(clin)$gene_name %in% c("TNFRSF10A", "TNFSF10"), ]
trail <- scale(t(trail), center = F)
trail <- rowSums(trail)
trail <- rank(trail)
trail <- trail/max(trail)
clin_df$trail <- trail
clin_df$trail.rank <- ifelse(trail < 0.30, "low", "hi")

# IFNG score, scale counts, sum counts, rank
ifng <- assay(clin)[rowRanges(clin)$gene_name %in% c("IFNG", "IFNGR1"), ]
ifng <- scale(t(ifng), center = F)
ifng <- rowSums(ifng)
ifng <- rank(ifng)
ifng <- ifng/max(ifng)
clin_df$ifng <- ifng
clin_df$ifng.rank <- ifelse(ifng < 0.30, "low", "hi")

png("Fig6_surv_dub_gbm.png", width = 1000, height = 750, res = 300)
ggsurvplot(survfit(Surv(overall_survival, deceased) ~ trail.rank == "hi" & ifng.rank == "hi", data = clin_df),
           palette = c("black", "green3"),
           censor.size = 3,
           font.title = c(10),
           font.x = c(7),
           font.y = c(7),
           font.tickslab = c(7),
           font.legend = c(7),
           legend = "none",
           xlab = c("Time (d)"),
           pval.size = 2.75,
           pval = T) +
  ggtitle("GBM")
dev.off()
