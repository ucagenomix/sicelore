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

custom.pal <- c("grey","#FF9B00","#EC380B")

markers <- c("Gad1","Gad2","Maf","Sst","Foxp1","Isl1","Htr3a","Dlx6os1","Dlx2","Nrxn3","Synpr","Slc17a6","Fabp7","Reln","Eomes","Vim","Top2a","Ube2c","Mki67","Ccnd2","Olig1","Pecam1","Igfbp7","Id2","Neurod6","Neurod2")
markers2 <- c("Olig1","Gad1","Gad2","Slc17a6","Slc17a7","Snap25","Gfap","Id2","Neurod1","Neurod2","Neurod6","Cd63","Rora","Calb2","Nr2f2","Htr3a","Dcx","Pax6","Eomes","Vim","S100b","Rbfox1","Ccnd2","Synpr","Mef2c","Pde1c")
precursors <- c("Mef2c","Erbb4","Plcxd3","Tspan7","Satb1","Synpr","Reln","Mpped1","Id2","Top2a","Cenpf","Ube2c","Olig1","Calb2","Pecam1","Cldn5","Igfbp7","Gad1","Gad2","Nrgn","Rora","Unc5c","Mdk","Neurod2","Slc17a6","Dlx1","Pbx3","Htr3a","Ckb","Cd63","Cd9","Slc6a5")
sub.markers <- c("Npy","Calb2","Nr2f2","Th","Cck","Reln","Pvalb","Nos1","Htr3a","Ccnd2","Id2","Synpr","Mef2c","Reln","Neurod6")

splicing.factor <- c('Celf1','Celf2','Celf4','Elavl4','Hnrnpa1','Hnrnpa2b1','Hnrnpa2b1','Hnrnpc','Hnrnpf','Hnrnph1','Hnrnph1','Hnrnph2','Hnrnpl','Khdrbs1','Khdrbs1','Khdrbs2','Khdrbs3','Mbnl1','Mbnl1','Mbnl2','Mbnl2','Mbnl3','Mbnl3','Nova1','Nova2','Ptbp1','Ptbp2','Ptbp3','Qk','Rbfox1','Rbfox2','Rbfox3','Rbm10','Rbm11','Rbm12','Rbm12b1','Rbm12b2','Rbm14','Rbm15','Rbm15b','Rbm17','Rbm18','Rbm19','Rbm20','Rbm22','Rbm24','Rbm25','Rbm26','Rbm27','Rbm28','Rbm3','Rbm31y','Rbm33','Rbm34','Rbm38','Rbm39','Rbm4','Rbm41','Rbm42','Rbm43','Rbm44','Rbm45','Rbm46','Rbm47','Rbm48','Rbm4b','Rbm5','Rbm6','Rbm7','Rbm8a','Sfpq','Srrm4','Srsf1','Srsf10','Srsf11','Srsf12','Srsf2','Srsf3','Srsf4','Srsf5','Srsf6','Srsf7','Srsf9','Tia1')

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 12)

setwd("/data/10x_data/10x_rainer/github")


