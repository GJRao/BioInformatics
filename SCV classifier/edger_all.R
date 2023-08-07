library(edgeR)
library(ggplot2)
library(ggrepel)


# get count dataset
counts <- read.csv("filtered_data_mrna.csv",row.names = 1L)
metadata <- read.csv("metadata.csv",row.names = 1L)
last_vec <- metadata[ , ncol(metadata)]                   # Apply ncol function
sample_info <- last_vec  
count_matrix <- as.matrix(counts)

dge <- DGEList(count = count_matrix, group = factor(sample_info))
dge
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

#calcNormFactors function is used for TMM normalization and calculating normalization (scaling) factors
dge <- calcNormFactors(object = dge)

#View dispersion estimates and biological coefficient of variation
jpeg("plotBCV.jpg")
plotBCV(dge)
dev.off()

dge <- estimateDisp(y = dge)

dge
#The exact test is performed using the exactTest() function and is applied only on a single factor design.
et <- exactTest(object = dge)
et
#topTags() function is useful to extract the table with adjusted p values (FDR). 
#The output table is ordered by p values.
#As the comparison of groups is Alive - Deceased, the positive log fold change represents the 
#gene is more highly expressed in the trt condition as compared to the ctr condition and vice versa.
top_degs <- topTags(object = et, n = Inf)
write.csv(top_degs,"top_degs.csv")
#Get a summary DGE table (returns significant genes with absolute log fold change at least 1 and adjusted p value < 0.05)
summary(decideTests(object = et,  adjust.method="BH", p.value=0.05, lfc=0.5))


degs_df <- as.data.frame(et$table)

#labelling genes
degs_df$diffexpressed <- "NO"
degs_df$diffexpressed[degs_df$logFC > 1 & degs_df$PValue < 0.05] <- "UP"
degs_df$diffexpressed[degs_df$logFC < -1 & degs_df$PValue < 0.05] <- "DOWN"

write.csv(degs_df,"degs_df1.csv")

# create an MA plot
#Make a MA plot of the libraries of count data
jpeg("plotMAResults1.jpg")
#plotSmear(et)
plotSmear(et,de.tags = rownames(et$table)[which(et$table$FDR>2)])
plotMA(et$table, bg.pch=1, bg.cex=0.5)
dev.off()

tab <- et$table
plot(tab$logCPM, tab$logFC, xlab="Average log CPM", ylab="log-fold-change", pch=16, cex=0.4)

#Make a MD plot of logFC against logcpm
jpeg("plotMDResults.jpg")
plotMD(et,bg.pch=1, bg.cex=0.5)
abline(h=c(-1.5, 1), col="blue")
dev.off()

