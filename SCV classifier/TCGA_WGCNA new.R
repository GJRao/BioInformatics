library(WGCNA)
library(DESeq2)
library(tidyverse)
library(cowplot)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

allowWGCNAThreads() 
#BiocManager::install("WGCNA")
data_exp <- read.csv("filtered_data_mrna.csv",row.names = 1L)
#phenodata
phenodata <- read.csv("metadata.csv",row.names = 1L)

traits <- read.csv("traits.csv",row.names = 1L)
#inverse
data_exp = t(data_exp)
#good sample genes
gsg <- goodSamplesGenes (data_exp, verbose=3)
gsg$allOK
summary(gsg)
table(gsg$goodGenes)
table(gsg$goodSamples)
#===============================================================================
#cluster samples to find outliers
sampleTree=hclust(dist(data_exp),method = "average")
#plot tree
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree)
rownames(data_exp)
pca <- prcomp(data_exp)
pca_data <- pca$x
#explained variance
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100,digits = 2)

pca_data <- as.data.frame(pca_data)
ggplot(pca_data,aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca_data))+
  labs(x=paste0('PC1:',pca.var.percent[1], ' %'),
       y=paste0('PC2:',pca.var.percent[2], ' %'))
#===============================================================================
#exclude outliers samples
sample.to.be.excluded <- c('TCGA.HZ.8003.01','TCGA.2J.AABV.01','TCGA.IB.8126.01',
                           'TCGA.LB.A9Q5.01','TCGA.IB.AAUM.01','TCGA.F2.6880.01')

data_exp <- t(data_exp)
data_exp.subset <- data_exp[,!(colnames(data_exp) %in% sample.to.be.excluded)]
#===============================================================================
#create Deseq2 dataset
#===============================================================================
colData <- phenodata %>%
  filter(!row.names(.) %in% sample.to.be.excluded)
names(colData) <- gsub('.S','_s',names(colData))

#colData <- traits
#making rownames and colnames identical
all(rownames(colData) %in% colnames(data_exp.subset))

#create dds
dds <- DESeqDataSetFromMatrix(countData = round(data_exp.subset),
                              colData = colData,
                              design = ~ Overall_survival_status)
#remove gene with <15 count
dds75 <- dds[rowSums(counts(dds)>=15) >=24,]  #16146 genes
#variance stabiliation
dds_norm <- vst(dds75)

#get normalized count
assay(dds_norm)
norm_counts <- assay(dds_norm) %>%
  t()
#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================
#choose soft-thresholding power
power <- c(c(1:10),seq(from=12, to=50,by=2))

#call the n/w topology analysis function
soft <- pickSoftThreshold(norm_counts, powerVector = power,
                          networkType = 'signed', verbose = 5)
soft.data <- soft$fitIndices

#visualize the power

a1 <- ggplot(soft.data,aes(Power,SFT.R.sq,label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  geom_hline(yintercept = 0.8,color='red')+
  
  labs(x='Power',y = 'scale free topology model fit, signed R^2')+
  theme_classic()

a2 <- ggplot(soft.data,aes(Power,mean.k.,label = Power))+
  geom_point()+
  geom_text(nudge_y = 0.1)+
  labs(x='Power',y = 'Mean connectivity')+
  theme_classic()

plot_grid(a1, a2, labels=c("A", "B"), ncol = 1, nrow = 2)

#convert matrix to numeric
norm_counts[] <- sapply(norm_counts,as.numeric)

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

datExpr <- norm_counts #normalized data
soft_power <- 10
# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = soft_power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

################################################################################
# Option 2b: 
TOM = TOMsimilarityFromExpr(datExpr, power = soft_power)
dissTOM = 1-TOM 
dim(dissTOM)

#===============================================================================
#
#  Construct modules (proceed with the genetree from option 2b)
#
#===============================================================================
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#dev.off()

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(9,9)
par(cex = 1.7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.3
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "5-merged_Module_Tree 0.4.pdf", width = 12, height = 9)  

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
write.csv(merge$oldMEs,file="oldMEs.csv");
write.csv(merge$newMEs,file="newMEs.csv")

#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================
# Heatmap of new module eigen-genes and sample trait (e.g. Zone)

#col_ann <- colData[,c(1)]
col_ann <- traits[,c(1,2)]
col_ann <- data.frame(col_ann)
row.names(col_ann)<-rownames(colData)
colnames(col_ann) <- c(colnames(colData))
#col_ann$Overall_survival_status <- as.factor(col_ann$Overall_survival_status)
#col_ann <- col_ann[order(col_ann$Overall_survival_status),]
#col_ann$sample_ID <- NULL
#head(col_ann)
ann_color <- list("col_ann" = c("Deceased" = "red",
                                "Alive" = "green"))
data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(data)
row.names(data)=colnames(data_exp.subset)
library(pheatmap)
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=10,
         annotation_row = col_ann, annotation_colors = ann_color)
#dev.off()

#=====================================================================================
#
#  PART 2: Correlation between gene modules and microbial traits (continuous data)
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

 # traits <- colData %>%
 #  mutate(survival_status=ifelse(grepl('Deceased',Overall_survival_status),0,1))%>%
 #  select(2)
 
#sorting
 traits <- traits[order(row.names(traits)), ]

 table(rownames(MEs) == rownames(traits))


# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
write.csv(moduleTraitCor,file="moduleTrait_correlation0.4_1.csv")
write.csv(moduleTraitPvalue,file="moduleTrait_pValue0.4_1.csv")

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-order0.4_1.pdf", width = 18, height = 15)
par(mar = c(10, 15, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = colnames(MEs),
               cex.lab = 1.5,
               cex.lab.x = 2,
               cex.lab.y = 1.8,
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =2,
               zlim = c(-0.5,0.5),
               main = paste("Module-trait relationships"),
               cex.main = 2.2)
dev.off()
#factors <-factor(colData$Overall_survival_status)
#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
#
#=====================================================================================


# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(traits$Deceased);
names(Verru) = "Deceased"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = (cbind(Verru,MEs))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
write.csv(MMPvalue,"Module membership p-val 0.4_1.csv")


geneTraitSignificance = as.data.frame(cor(datExpr, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="")
write.csv(GSPvalue,"Gene significance p-val0.4_1.csv")
# Plot the dendrogram
sizeGrWindow(6,6)
pdf("eigengenedendo.pdf", width = 15, height = 10)
par(mar = c(10, 15, 5, 5));
par(cex = 2.5)
plotEigengeneNetworks(MET, "Eigengene dendrogram1", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf("6-Eigengene adjacency heatmap 0.4_1.pdf", width = 20, height = 15)
par(cex = 2.5)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(8,8,1,1),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

module = "purple"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf("scatter plot.pdf", width = 15, height = 10)
sizeGrWindow(6,6);

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Deceased",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 2, cex.lab = 1.5, cex.axis = 1.5, col = "blue")

dev.off()
#module color mapping to gene
mod_gen_mapping <- cbind(rownames(geneModuleMembership), as.data.frame(merge$colors))
colnames(mod_gen_mapping) <- c("Genes","Module")
write.csv(mod_gen_mapping,"Module-gene-mapping0.4_1.csv")


blue_mod_genes <- mod_gen_mapping[mod_gen_mapping$Module=="blue",] 
write.csv(blue_mod_genes,"blue_mod_genes.csv")
# ## Draw bubble plot for particular module
# colsum_traits <- colSums(traits)
# colsum_traits <- data.frame(colsum_traits)
# colsum_traits$b_order <- rownames(colsum_traits)
# library(tidyr)
# moduleTraitCor_long <- data.frame(moduleTraitCor)
# moduleTraitCor_long$module <- rownames(moduleTraitCor)
# moduleTraitCor_long <- moduleTraitCor_long[,c(2)]
# moduleTraitCor_long <- gather(moduleTraitCor_long, b_order, PCC, Pseudomonadales:Others, factor_key = TRUE)

