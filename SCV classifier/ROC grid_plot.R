
library(png)
library(raster)
library(tidyverse)
library(ggplot2)
library(cowplot)

library(gridExtra)
library(ggpubr)

img0 <- readPNG("AURKB1.png")
img0_grob <- grid::rasterGrob(img0)

img1 <- readPNG("bub1_1.png")
img1_grob <- grid::rasterGrob(img1)

img2 <- readPNG("bub1b_1.png")
img2_grob <- grid::rasterGrob(img2)

img3 <- readPNG("ccna2_1.png")
img3_grob <- grid::rasterGrob(img3)

img4 <- readPNG("kif11_1.png")
img4_grob <- grid::rasterGrob(img4)


img5 <- readPNG("KIF23_1.png")
img5_grob <- grid::rasterGrob(img5)

img6 <- readPNG("TOP2A_1.png")
img6_grob <- grid::rasterGrob(img6)

img7 <- readPNG("TTK_1.png")
img7_grob <- grid::rasterGrob(img7)

roc_grid <- plot_grid(img0_grob,img1_grob,img2_grob, 
                                 img3_grob,img4_grob, img5_grob,
                                 img6_grob,img7_grob,
                                 labels = c("A","B","C","D","E","F","G","H"),
                                 ncol=3)
                      
ggsave("roc_grid3.jpeg",roc_grid,dpi=400,width = 10, height = 11)

plot_grid(img0_grob,img1_grob)


