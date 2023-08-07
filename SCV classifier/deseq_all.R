library(DESeq2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggfortify)

# get count dataset
counts <- read.csv("filtered_data_mrna.csv",row.names = 1L)
metadata <- read.csv("metadata.csv",row.names = 1L)

cols <- colnames(counts)
#check row names in meta equals to cols names in counts data
all(cols %in% rownames(metadata))
#colnames of counts_data and row names of meta data are in same order?
all(cols == rownames(metadata))

#create DESEQ obj
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Overall.Survival.Status)
dds
#removing genes with low count.keeping genes with atleast 10 row count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#set the factor level. normal vs tumor in metadata
dds$Condition <- relevel(dds$Overall.Survival.Status, ref = "Deceased")
dds <- DESeq(dds)

# see all comparisons (here there is only one)
resultsNames(dds)
# get gene expression table
# at this step independent filtering is applied by default to remove low count genes
# independent filtering can be turned off by passing independentFiltering=FALSE to results
res <- results(dds)  # same as results(dds, name="condition_infected_vs_control") or results(dds, contrast = c("condition", "infected", "control") )
res
#Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method) 
res[order(res$padj),] 

write.csv(as.data.frame(res[order(res$padj),] ), file="Overall.Survival.Status Deceased vs Alive_dge.csv")
#Get summary of differential gene expression with adjusted p value cut-off at 0.05,
res<-results(dds, alpha=0.01)
summary(res)

#Normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)

#BiocManager::install("apeglm")

#Shrinkage estimation of log2 fold changes (LFCs)
resLFC <- lfcShrink(dds, coef="Overall.Survival.Status_Deceased_vs_Alive", type="apeglm")
head(resLFC)

par(mfrow = c(1, 2))
plotMA(resLFC, main="Shrinkage of LFCs", ylim=c(-4,4))
plotMA(res, main="No shrinkage of LFCs", ylim=c(-4,4))

#MAplot
par(mfrow = c(1, 1))
plotMA(res)


#labelling genes
res_frame <- as.data.frame(res)
res_frame$diffexpressed <- "NO"
res_frame$diffexpressed[res_frame$log2FoldChange > 1 & res_frame$pvalue < 0.05] <- "UP"
res_frame$diffexpressed[res_frame$log2FoldChange < -1 & res_frame$pvalue < 0.05] <- "DOWN"
write.csv(res_frame,'DEGs_with_updown_pval1.csv')
