rm(list = ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)

scRNA <- readRDS('data/scRNA.rds')
Idents(scRNA) <- 'celltype'
table(scRNA$celltype)

#T/nk细胞
T_cell <- subset(scRNA,idents = 'T/NK')

ElbowPlot(T_cell, ndims=50, reduction="harmony")
pct <- T_cell[["pca"]]@stdev / sum(T_cell[["harmony"]]@stdev) * 100 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use

T_cell <- FindNeighbors(T_cell, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = seq(0.4,1,by=0.1))
T_cell <- RunTSNE(T_cell, reduction = "harmony", dims = 1:20)
T_cell <- RunUMAP(T_cell, reduction = "harmony", dims = 1:20)
table(T_cell$RNA_snn_res.0.5)
Idents(T_cell) <- 'RNA_snn_res.0.5'


gene <- c('CD3E','CD8A','CD4','NCAM1','CD34','CCR7','IL7R','SELL','TCF7','CD44','TIGIT','HAVCR2','CTLA4','LAG3','CD274',
          'ICOS','ICOSL','CD28','TNFRSF18','TNFRSF8','CD40LG','FOXP3','IKZF2','PDCD1','CXCR5','CXCL13','ADGRG1','MAF','PRDM1',
          'SOX4','BCL6','BATF','KLRB1','CCR4','TBX21','GATA3','IL2','IL4','IFNG','TRP53','EGFR','PTEN','MYC','DNMT3A','TET2','RHOA',
          'MS4A1','CD19','TNFRSF1B','UBA52','SATB1','PRF1','NKG7','GZMA','GZMB','KLRG1','ITGAE','CIITA','TGFB1','CD38','NFKB2','AR','CCR2','CCR5',
          'IL21','CX3CR1','RORA','STAT3','CCR6','IL1R1','IL23R','IL17F','MKI67','MME','TTBK1','HSPA1A','RORB','NCR1','IL2RA','CXCR3','CD24A','SDC1',
          'CD27','CX3CR1','IL17A','CR1','CR2')
DotPlot(T_cell, features = unique(gene),
        assay='RNA'  )  + coord_flip()+
  scale_color_gradientn(colours = viridis::viridis(20))
ggsave(filename = 'celltype/Tcell_10cluster.pdf',width = 10,height = 15)

saveRDS(T_cell,'data/AITL_ONLY_T_NK.rds')
table(T_cell$celltype)

T_cell <- readRDS('data/AITL_ONLY_T_NK.rds')
DimPlot(T_cell,reduction = 'harmony',group.by = 'orig.ident')
DimPlot(scRNA,reduction = 'harmony',group.by = 'orig.ident')
#B
scRNA <- readRDS('data/scRNA.rds')
Idents(scRNA) <- 'celltype'
table(scRNA$celltype)
B_cell <- subset(scRNA,idents = 'B')

ElbowPlot(B_cell, ndims=50, reduction="harmony")
pct <- B_cell[["pca"]]@stdev / sum(B_cell[["harmony"]]@stdev) * 100 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use

B_cell <- FindNeighbors(B_cell, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = seq(0.4,1,by=0.1))
B_cell <- RunTSNE(B_cell, reduction = "harmony", dims = 1:20)
B_cell <- RunUMAP(B_cell, reduction = "harmony", dims = 1:20)
table(B_cell$RNA_snn_res.0.5)
table(B_cell$celltype)
saveRDS(B_cell,file = 'data/AITL_ONLT_B.rds')

#Myeloid
Myeloid <- subset(scRNA,idents = 'Myeloid')

ElbowPlot(Myeloid, ndims=50, reduction="harmony")
pct <- Myeloid[["pca"]]@stdev / sum(Myeloid[["harmony"]]@stdev) * 100 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use

Myeloid <- FindNeighbors(Myeloid, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = seq(0.4,1,by=0.1))
Myeloid <- RunTSNE(Myeloid, reduction = "harmony", dims = 1:20)
Myeloid <- RunUMAP(Myeloid, reduction = "harmony", dims = 1:20)
table(Myeloid$RNA_snn_res.0.5)
table(Myeloid$celltype)
saveRDS(Myeloid,file = 'data/AITL_ONLY_Myeloid.rds')
#stromal
stromal <- subset(scRNA,idents = 'stromal')

ElbowPlot(stromal, ndims=50, reduction="harmony")
pct <- stromal[["pca"]]@stdev / sum(stromal[["harmony"]]@stdev) * 100 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use

stromal <- FindNeighbors(stromal, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = seq(0.4,1,by=0.1))
stromal <- RunTSNE(stromal, reduction = "harmony", dims = 1:20)
stromal <- RunUMAP(stromal, reduction = "harmony", dims = 1:20)
table(stromal$RNA_snn_res.0.5)
table(stromal$celltype)
saveRDS(stromal,file = 'data/AITL_ONLY_stromal.rds')

#Check
list.files('data')
T_cell <- readRDS('data/AITL_ONLY_T_NK.rds')
B_cell <- readRDS('data/AITL_ONLY_B.rds')
Myeloid <- readRDS('data/AITL_ONLY_Myeloid.rds')
stromal <- readRDS('data/AITL_ONLY_stromal.rds')
table(T_cell$celltype)
table(B_cell$celltype)
table(Myeloid$celltype)
table(stromal$celltype)
