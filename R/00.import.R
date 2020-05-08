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
library(DropletUtils)
library(reshape2)
library(gameofthrones)
library(data.table)
library(plyr)

custom.pal <- c("grey","#FF9B00","#EC380B")
cell_type_color <- c("radial glia"="#E54532","cycling radial glia"="#7A5388","intermediate progenitor"="#2A76BB","imature Glutamatergic"="#38AEEF","mature Glutamatergic"="#4CB57C","imature GABAergic"="#A5A739","mature GABAergic"="#EA766B","Cajal-Retzius"="#F0E816")
clusters <- c("radial glia","cycling radial glia","intermediate progenitor", "imature Glutamatergic", "mature Glutamatergic","imature GABAergic","mature GABAergic", "Cajal-Retzius")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
col=ggplotColours(n = 7)


# set scale.data slot from data slot
scale_my_data <- function(seurat, assay){
  data_in <- as.matrix(seurat@assays[[assay]]@data)
  mat <- sapply(strsplit(rownames(data_in), "\\.\\."), `[`, 1)
  
  all.genes <- unique(mat)
  all.genes <- as.data.frame(all.genes)
  
  print(paste("Unique features =", dim(all.genes)[1], "out of", dim(data_in)[1], "total features", sep=" "))
  
  colnames(all.genes)<-c("geneId")
  rownames(all.genes) <- all.genes$geneId
  
  output <- data_in
  
  # for all genes
  for(i in 1:dim(all.genes)[1]){
    x <- which(mat == all.genes[i,])
    m <- mean(unlist(data_in[x,]))
    sd <- sd(unlist(data_in[x,]))
    
    # for all gene isoforms
    for(j in 1:length(x)){
      
      # for all cells
      for(k in 1:dim(data_in)[2]){
        output[x[j],k] <- (data_in[x[j],k] - m)/sd
      }
    }
  }
  
  seurat <- SetAssayData(seurat, slot = "scale.data", output, assay = assay)
  return(seurat)
}

                  