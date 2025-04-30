rm(list = ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
scRNA <- readRDS('data/scRNA.rds')
Idents(scRNA) <- 'RNA_snn_res.0.5'

table(scRNA@active.ident)

gene <- c('CD2','CD3D','CD3G','TRAC','IL32',#T/NK
          'IGFBP7','SPARCL1','MGP','COL1A2','DCN',#stromal
          'LYZ','AIF1','TYROBP','FCER1G','IGSF6',#Myeloid
          'MS4A1','CD79A','IGHM','IGHG3','IGHG1'#B cells
          )
DotPlot(scRNA, features = unique(gene),
        assay='RNA'  )  + coord_flip()+
  scale_color_gradientn(colours = viridis::viridis(20))
ggsave(filename = 'celltype/Celltype_annotation.pdf',width = 8,height = 6)

celltype=data.frame(ClusterID=0:15,
                    celltype= 0:15)
celltype[celltype$ClusterID %in% c(0,1,2,3,5,8,10,11,12),2]='T/NK'
celltype[celltype$ClusterID %in% c(15),2]='stromal'
celltype[celltype$ClusterID %in% c(6,13,14),2]='Myeloid'
celltype[celltype$ClusterID %in% c(4,7,9),2]='B'
table(celltype$celltype)

scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype,scRNA$RNA_snn_res.0.5)

saveRDS(scRNA,file = 'data/scRNA.rds')
