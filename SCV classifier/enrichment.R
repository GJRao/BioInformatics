library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

 gene_list <- read.csv("SCV_res_deg_updown.csv")
 purplelist <- read.csv("scv_wgcna_purple.csv")
# # 
final.gene.list <- gene_list[gene_list$genes %in% purplelist$genes,]
write.csv(final.gene.list,"scv_wgcna_purple with up down.csv")
# write.csv(final.gene.list,"scv_wgcna_purple.csv")
BiocManager::install("clusterProfiler")

library(DOSE)
#BiocManager::install("org.Hs.eg.db")


d <- read.csv("scv_wgcna_purple.csv")

d$ENTREZID <- mapIds(org.Hs.eg.db, d$genes, 'ENTREZID', 'SYMBOL')

#Remove rows with NA's using na.omit()
d <- na.omit(d)

#feature 1:foldchange column 4
genelist <- d[,4]


## feature 2: named vector
names(genelist) = as.character(d[,9])


## feature 3: decreasing orde
genelist = sort(genelist, decreasing = TRUE)


#de <- names(genelist)[abs(genelist) < 0.01]
de <- names(genelist)[1:100]


#Gene Ontology enrichment

ego1 <- enrichGO(gene          = de,
                universe      = names(genelist),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
dotplot(ego1, showCategory=10)+ ggtitle("GO_BP")

ego2 <- enrichGO(gene          = de,
                 universe      = names(genelist),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                readable      = TRUE)
head(ego)
dotplot(ego2, showCategory=10)+ ggtitle("GO_CC")

ego3 <- enrichGO(gene          = de,
                 universe      = names(genelist),
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                readable      = TRUE)
head(ego)
dotplot(ego3, showCategory=10)+ ggtitle("GO_MF")

goplot(ego)

#KEGG pathway
#KEGG pathway over-representation analysis

kk <- enrichKEGG(gene         = de,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

dotplot(kk, showCategory=10)+ ggtitle("KEGG pathway")
goplot(kk)

#Gene-Concept Network

x <- enrichDO(de)
cnetplot(x)

#Gene-Concept Network
edo <- enrichDGN(de)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=genelist)

cnetplot(edox, foldChange=genelist, circular = TRUE, colorEdge = TRUE)



## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=genelist)