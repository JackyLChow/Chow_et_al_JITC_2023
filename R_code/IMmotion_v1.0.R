library(ggplot2)
library(gridExtra)
library(dplyr)
library(readr)
library(stringr)

# counts data EGAD00001006618, log2(TPM + 1)
ct <- read_csv("Data/IMmotion/IMmotion151.expression.data.TPM.anon.20201106.csv")
# clinical metadata EGAD00001006617
cd <- read_csv("Data/IMmotion/IMmotion151_clinical_anon_20201106_trim.csv")

# filter meta data for progression evaluated tumors and available RNAseq data
meta <- cd[cd$OBJECTIVE_RESPONSE != "NE" &
             !is.na(cd$RNASEQ_SAMPLE_ID), ]
# filter counts columns
counts <- ct[, c("symbol", meta$RNASEQ_SAMPLE_ID)]

# bind data
meta$RNASEQ_SAMPLE_ID <- str_replace(meta$RNASEQ_SAMPLE_ID, "-", ".")

genes <- c("IFNG", "IFNGR1", "TNFSF10", "TNFRSF10A")
genes <- data.frame(counts[counts$symbol %in% genes, ])
rownames(genes) <- genes$symbol
genes <- t(as.matrix(genes[, 2:ncol(genes)]))

data <- cbind(meta, genes)

# remove sarcomatoid, primary tumors only, sort by arm
## atezo_bev
foo <- filter(data, ARM == "atezo_bev", SARCOMATOID == "no", PRIMARY_VS_METASTATIC == "PRIMARY")
foo$CR <- ifelse(foo$OBJECTIVE_RESPONSE == "CR", T, F) # separate CR from others
foo$OBJECTIVE_RESPONSE <- factor(foo$OBJECTIVE_RESPONSE, levels = c("CR", "PR", "SD", "PD")) # establish levels

anova(lm(foo$TNFSF10 + foo$TNFRSF10A + foo$IFNG + foo$IFNGR1 ~ foo$OBJECTIVE_RESPONSE)) # anova stat test

png("Fig6_atez.png", width = 1000, height = 750, res = 300)
ggplot(foo, aes(x = OBJECTIVE_RESPONSE, y = (IFNG + IFNGR1 + TNFSF10 + TNFRSF10A) * 4 / 100)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(title = "atezo_bev arm\np = 0.01",
       y = "Pre-treatment IFNG/R TRAIL/R") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
dev.off()

t.test(foo$TNFSF10 + foo$TNFRSF10A + foo$IFNG + foo$IFNGR1 ~ foo$CR)
t.test(foo$IFNG ~ foo$CR)
t.test(foo$IFNGR1 ~ foo$CR)
t.test(foo$TNFSF10 ~ foo$CR)
t.test(foo$TNFRSF10A ~ foo$CR)

png("SFig7_cr_score.png", width = 500, height = 750, res = 300)
ggplot(foo, aes(x = CR, y = (IFNG + IFNGR1 + TNFSF10 + TNFRSF10A) * 4 / 100)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_x_discrete(limits = c(T, F),
                   labels = c("CR", "PR, SD, PD")) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(title = "p = 0.002",
       y = "Pre-treatment IFNG/R TRAIL/R") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
dev.off()

png("SFig7_cr_ind.png", width = 1750, height = 750, res = 300)
grid.arrange(grobs = list(
ggplot(foo, aes(x = CR, y = IFNG)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_x_discrete(limits = c(T, F),
                   labels = c("CR", "PR, SD, PD")) +
  labs(title = "p = 0.14",
       y = "IFNG normalized counts") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
,
ggplot(foo, aes(x = CR, y = IFNGR1)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_x_discrete(limits = c(T, F),
                   labels = c("CR", "PR, SD, PD")) +
  labs(title = "p = 0.68",
       y = "IFNGR1 normalized counts") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
,
ggplot(foo, aes(x = CR, y = TNFSF10)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_x_discrete(limits = c(T, F),
                   labels = c("CR", "PR, SD, PD")) +
  labs(title = "p = 0.006",
       y = "TNFSF10 normalized counts") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
,
ggplot(foo, aes(x = CR, y = TNFRSF10A)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_x_discrete(limits = c(T, F),
                   labels = c("CR", "PR, SD, PD")) +
  labs(title = "p = 0.10",
       y = "TNFRSF10A normalized counts") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "black"))
),
nrow = 1)
dev.off()

## sunitinib
foo <- filter(data, ARM == "sunitinib", SARCOMATOID == "no", PRIMARY_VS_METASTATIC == "PRIMARY")
foo$OBJECTIVE_RESPONSE <- factor(foo$OBJECTIVE_RESPONSE, levels = c("CR", "PR", "SD", "PD"))

anova(lm(foo$TNFSF10 + foo$TNFRSF10A + foo$IFNG + foo$IFNGR1 ~ foo$OBJECTIVE_RESPONSE))

png("Fig6_sunit.png", width = 500, height = 750, res = 300)
ggplot(foo, aes(x = OBJECTIVE_RESPONSE, y = (IFNG + IFNGR1 + TNFSF10 + TNFRSF10A) * 4 / 100)) +
  geom_boxplot(color = "black", 
               outlier.size = 0.5) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(title = "sunitinib arm\np = 0.11",
       y = "Pre-treatment IFNG, TRAIL score") +
  stat_summary(fun = "mean", color = "black", size = 0.5, shape = 23, fill = "white") +
  theme_classic() +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y = element_text(size = 7, color = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7, color = "white"))
dev.off()

