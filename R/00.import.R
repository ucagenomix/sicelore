library(Seurat)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(knitr)
library(xtable)
library(ggrepel)
library(cowplot)
library(metap)
library(DoubletFinder)

source("/data/10x_data/DoubletFinder-master/R/paramSweep_v3.R")
source("/data/10x_data/DoubletFinder-master/R/doubletFinder_v3.R")
source("/data/10x_data/DoubletFinder-master/R/summarizeSweep.R")
source("/data/10x_data/DoubletFinder-master/R/find.pK.R")
source("/data/10x_data/DoubletFinder-master/R/modelHomotypic.R")

custom.pal <- c("grey","#FF9B00","#EC380B")

#cc.genes <- readLines(con = "/data/10x_data/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_mouse.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
#s.genes <- cc.genes[1:43]
#g2m.genes <- cc.genes[44:97]


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 12)
