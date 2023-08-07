
library(png)
library(raster)
library(tidyverse)
library(ggplot2)
library(cowplot)

library(gridExtra)
library(ggpubr)

img0 <- readPNG("AURKB.png")
img0_grob <- grid::rasterGrob(img0)

img1 <- readPNG("BUB1.png")
img1_grob <- grid::rasterGrob(img1)

img2 <- readPNG("BUB1B.png")
img2_grob <- grid::rasterGrob(img2)

img3 <- readPNG("CCNA2.png")
img3_grob <- grid::rasterGrob(img3)

img4 <- readPNG("KIF11.png")
img4_grob <- grid::rasterGrob(img4)


img5 <- readPNG("KIF23.png")
img5_grob <- grid::rasterGrob(img5)

img6 <- readPNG("TOP2A.png")
img6_grob <- grid::rasterGrob(img6)

img7 <- readPNG("TTK.png")
img7_grob <- grid::rasterGrob(img7)

gepia_grid <- plot_grid(img0_grob,img1_grob,img2_grob, 
                                 img3_grob,img4_grob, img5_grob,
                                 img6_grob,img7_grob,
                                 labels = c("A","B","C","D","E","F","G","H"),
                                 ncol=4, nrow = 2)
                      
ggsave("gepia_grid.jpeg",gepia_grid,dpi =700,width = 10, height = 8)




