library(tidyverse)
library(pROC)
library(cowplot)

data_exp <- read.csv("filtered_data_mrna.csv",row.names = 1L)

hub_gene <- c("AURKB", "BUB1", "BUB1B", "CCNA2", "KIF11", "KIF23",
              "TOP2A","TTK")  

#hub_gene1 <- c("EXO1", "AURKA", "CCNB1", "CDC20")
data_exp1 <- data_exp[rownames(data_exp) %in% hub_gene1,]

data_exp_trans <- t(data_exp1)
metadata <- read.csv("metadata.csv",row.names = 1L)
#data(data_exp)
surv <- metadata$Survival_status
data_exp_trans <- cbind(data_exp_trans, surv)
data_exp_trans <- as.data.frame(data_exp_trans)


#ROC curve object 

roc_AURKB <- roc(surv ~ AURKB, data_exp_trans,
                plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
                main="AURKB",cex.main = 1.5, cex.text=2, 
                cex.axis = 1.3,smooth = TRUE)

plot(roc_AURKB, cex.lab = 1.5, 
     main="AURKB (AUC= 0.624)",cex.main = 1.5, cex.text=2, cex.axis = 1.3,
     print.auc=TRUE, col="black",  lwd=2, 
     print.auc.y=0.45, cex.text= 2)



roc_BUB1 <- roc(surv ~ BUB1, data_exp_trans,
               plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
               main="BUB1 (AUC= 0.673)",cex.main = 1.5, cex.text=2, 
               cex.axis = 1.3,smooth = TRUE)

roc_BUB1B <- roc(surv ~ BUB1B, data_exp_trans,
               plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
               main="BUB1B (AUC= 0.663)",cex.main = 1.5, cex.text=2, 
               cex.axis = 1.3,smooth = TRUE)

roc_CCNA2 <- roc(surv ~ CCNA2, data_exp_trans,
               plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
               main="CCNA2 (AUC= 0.643)",cex.main = 1.5, cex.text=2, 
               cex.axis = 1.3,smooth = TRUE)

roc_KIF11 <- roc(surv ~ KIF11, data_exp_trans,
               plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
               main="KIF11 (AUC= 0.665)",cex.main = 1.5, cex.text=2, 
               cex.axis = 1.3,smooth = TRUE)

roc_KIF23 <- roc(surv ~ KIF23, data_exp_trans,
               plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
               main="KIF23 (AUC= 0.688)",cex.main = 1.5, cex.text=2, 
               cex.axis = 1.3,smooth = TRUE)



roc_TOP2A <- roc(surv ~ TOP2A, data_exp_trans,
                 plot = TRUE, print.auc = TRUE,cex.lab = 1.5, 
                 main="TOP2A (AUC= 0.682)",cex.main = 1.5, cex.text=2, 
                 cex.axis = 1.3,smooth = TRUE)


  roc_TTK <- roc(data_exp_trans$surv, data_exp_trans$TTK,
                  plot = TRUE, print.auc = TRUE, cex.lab = 1.5, 
                 main="TTK (AUC= 0.663)",cex.main = 1.5, cex.text=2, 
                 cex.axis = 1.3,smooth = TRUE)


