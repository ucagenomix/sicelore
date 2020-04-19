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

#cc.genes <- readLines(con = "/data/10x_data/cell_cycle_vignette_files/regev_lab_cell_cycle_genes_mouse.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
#s.genes <- cc.genes[1:43]
#g2m.genes <- cc.genes[44:97]

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
    
    #print(paste("i",i,sep="="))
    x <- which(mat == all.genes[i,])
    m <- mean(unlist(data_in[x,]))
    sd <- sd(unlist(data_in[x,]))
    
    # for all gene isoforms
    for(j in 1:length(x)){
      #print(paste("j",j,length(x),sep="="))
      
      # fo all cells
      for(k in 1:dim(data_in)[2]){
        #print(paste("k",k,sep="="))
        output[x[j],k] <- (data_in[x[j],k] - m)/sd
      }
    }
  }
  
  #seurat@assays[[assay]]@scale.data <- output
  seurat <- SetAssayData(seurat, slot = "scale.data", output, assay = assay)
  return(seurat)
}

                  