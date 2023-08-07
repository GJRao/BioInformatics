
library(png)
library(raster)
library(tidyverse)
library(ggplot2)
library(cowplot)

library(gridExtra)
library(ggpubr)

img0 <- readPNG("GO_BP.png")
img0_grob <- grid::rasterGrob(img0)

img1 <- readPNG("GO_CC.png")
img1_grob <- grid::rasterGrob(img1)

img2 <- readPNG("GO_MF.png")
img2_grob <- grid::rasterGrob(img2)

img3 <- readPNG("kegg.png")
img3_grob <- grid::rasterGrob(img3)

img4 <- readPNG("purple_SCV venn.png")
img4_grob <- grid::rasterGrob(img4)

img5 <- readPNG("enrichDGN.png")
img5_grob <- grid::rasterGrob(img5)


GO_grid <- plot_grid(img4_grob,img0_grob,img1_grob,img2_grob, img3_grob,img5_grob,
                       labels = c("A","B","C","D","E","F"),ncol=2, nrow = 3)
ggsave("GO_grid1.jpeg",GO_grid,dpi =700,width = 10, height = 10)




