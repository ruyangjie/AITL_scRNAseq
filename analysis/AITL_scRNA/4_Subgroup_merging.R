rm(list = ls())
gc()

scRNA <- readRDS('data/scRNA.rds')
T_cell <- readRDS('data/AITL_ONLY_T_NK.rds')
B_cell <- readRDS('data/AITL_ONLY_B.rds')

metadata <- scRNA@meta.data
T_data <- as.data.frame(T_cell$RNA_snn_res.0.5)
metadata$cellchat <- 'NA'
match(rownames(T_data),rownames(metadata)) 
for(i in 1:nrow(T_data)){
  metadata[which(rownames(metadata) == rownames(T_data)[i]),'cellchat'] <- paste0('T_',T_data$`T_cell$RNA_snn_res.0.5`[i])}
#check
identical(rownames(T_data),rownames(metadata[metadata$cellchat!= 'NA',]))
table(T_cell$RNA_snn_res.0.5)
table(metadata$cellchat)



B_data <- as.data.frame(B_cell$RNA_snn_res.0.4)
table(B_data$`B_cell$RNA_snn_res.0.4`)
match(rownames(B_data),rownames(metadata)) 
for(i in 1:nrow(B_data)){
  metadata[which(rownames(metadata) == rownames(B_data)[i]),'cellchat'] <- paste0('B_',B_data$`B_cell$RNA_snn_res.0.4`[i])}
table(metadata$cellchat)
table(B_cell$RNA_snn_res.0.4)
table(T_cell$RNA_snn_res.0.5)

scRNA@meta.data <- metadata

DimPlot(scRNA,reduction = 'umap',group.by = 'celltype')
Idents(scRNA) <- 'cellchat'
Tcell <- c('T_0','T_1','T_2','T_3','T_4','T_5','T_6','T_7','T_8','T_9')
DimPlot(scRNA,reduction = 'umap',group.by = 'cellchat',cols = c('T_0'='red','T_1'='red','T_2'='red','T_3'='red','T_4'='red','T_5'='red','T_6'='red','T_7'='red','T_8'='red','T_9'='red'))
DimPlot(scRNA,reduction = 'umap',group.by = 'cellchat',cols = c('B_0'='red','B_1'='red','B_2'='red','B_3'='red','B_4'='red','B_5'='red','B_6'='red'))

saveRDS(scRNA,file = 'data/cellchat_use_data.rds')

####################################
#annotation
####################################
rm(list = ls())
gc()
scRNA <- readRDS('data/scRNA.rds')
metadata <- scRNA@meta.data

####################################
#T
####################################
T_cell <- readRDS('data/AITL_ONLY_T_NK.rds')
T_data <- as.data.frame(T_cell$RNA_snn_res.0.5)
metadata$sub_celltype <- 'NA'
match(rownames(T_data),rownames(metadata)) 
table(T_cell$RNA_snn_res.0.5)
celltype=data.frame(ClusterID=0:9,
                    celltype= 0:9)
celltype[celltype$ClusterID %in% c(0),2]='γδ_T'
celltype[celltype$ClusterID %in% c(1),2]='TFH1'
celltype[celltype$ClusterID %in% c(2),2]='TFH2'
celltype[celltype$ClusterID %in% c(3),2]='DNT'
celltype[celltype$ClusterID %in% c(4),2]='TPH'
celltype[celltype$ClusterID %in% c(5,6),2]='TEMRA_CD8_T'
celltype[celltype$ClusterID %in% c(7),2]='Treg'
celltype[celltype$ClusterID %in% c(8),2]='Ki67+_γδ_T'
celltype[celltype$ClusterID %in% c(9),2]='TFH3'
table(celltype$celltype)
T_data$celltype = "NA"
for(i in 1:nrow(celltype)){
  T_data[which(T_data$`T_cell$RNA_snn_res.0.5` == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(T_data$celltype,T_data$`T_cell$RNA_snn_res.0.5`)

####################################
#B
####################################
B_cell <- readRDS('data/AITL_ONLY_B.rds')
B_data <- as.data.frame(B_cell$RNA_snn_res.0.4)
match(rownames(B_data),rownames(metadata))
table(B_cell$RNA_snn_res.0.4)

celltype=data.frame(ClusterID=0:6,
                    celltype= 0:6)
celltype[celltype$ClusterID %in% c(0),2]='naïve B'
celltype[celltype$ClusterID %in% c(1,4),2]='Plasma'
celltype[celltype$ClusterID %in% c(2,3),2]='Plasmablast'
celltype[celltype$ClusterID %in% c(5),2]='DPT'
celltype[celltype$ClusterID %in% c(6),2]='T_B'
table(celltype$celltype)
B_data$celltype = "NA"
for(i in 1:nrow(celltype)){
  B_data[which(B_data$`B_cell$RNA_snn_res.0.4` == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(B_data$celltype,B_data$`B_cell$RNA_snn_res.0.4`)

####################################
#myeloid
####################################
M_cell <- readRDS('data/AITL_ONLY_Myeloid.rds')
M_data <- as.data.frame(M_cell$RNA_snn_res.0.4)
match(rownames(M_data),rownames(metadata))
table(M_cell$RNA_snn_res.0.4)

celltype=data.frame(ClusterID=0:6,
                    celltype= 0:6)
celltype[celltype$ClusterID %in% c(0),2]='pDC'
celltype[celltype$ClusterID %in% c(1),2]='Monocytes'
celltype[celltype$ClusterID %in% c(2),2]='Neutrophils'
celltype[celltype$ClusterID %in% c(3),2]='Macrophage(M2)'
celltype[celltype$ClusterID %in% c(4),2]='TFH3'
celltype[celltype$ClusterID %in% c(5),2]='undefined'
celltype[celltype$ClusterID %in% c(6),2]='cDC'
table(celltype$celltype)
M_data$celltype = "NA"
for(i in 1:nrow(celltype)){
  M_data[which(M_data$`M_cell$RNA_snn_res.0.4` == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(M_data$celltype,M_data$`M_cell$RNA_snn_res.0.4`)

####################################
#stromal
####################################
S_cell <- readRDS('data/AITL_ONLY_stromal.rds')
S_data <- as.data.frame(S_cell$RNA_snn_res.0.4)
match(rownames(S_data),rownames(metadata))
table(S_cell$RNA_snn_res.0.4)
celltype=data.frame(ClusterID=0:4,
                    celltype= 0:4)
celltype[celltype$ClusterID %in% c(0),2]='FDC'
celltype[celltype$ClusterID %in% c(1,3),2]='endothelial_cells'
celltype[celltype$ClusterID %in% c(2,4),2]='T_FDC'
table(celltype$celltype)
S_data$celltype = "NA"
for(i in 1:nrow(celltype)){
  S_data[which(S_data$`S_cell$RNA_snn_res.0.4` == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(S_data$celltype,S_data$`S_cell$RNA_snn_res.0.4`)

#Merge
####################################
#T
####################################
for(i in 1:nrow(T_data)){
  metadata[which(rownames(metadata) == rownames(T_data)[i]),'sub_celltype'] <- T_data$celltype[i]}

identical(rownames(T_data),rownames(metadata[metadata$sub_celltype!= 'NA',]))#顺序和内容都一样
####################################
#B
####################################
for(i in 1:nrow(B_data)){
  metadata[which(rownames(metadata) == rownames(B_data)[i]),'sub_celltype'] <- B_data$celltype[i]}

identical(rownames(B_data),rownames(metadata[metadata$sub_celltype%in%B_data$celltype,]))#顺序和内容都一样

####################################
#M
####################################
for(i in 1:nrow(M_data)){
  metadata[which(rownames(metadata) == rownames(M_data)[i]),'sub_celltype'] <- M_data$celltype[i]}

identical(rownames(M_data),rownames(metadata[metadata$sub_celltype%in%M_data$celltype,]))#顺序和内容都一样
#这里混进去了T的细胞所以不一样

####################################
#S
####################################
for(i in 1:nrow(S_data)){
  metadata[which(rownames(metadata) == rownames(S_data)[i]),'sub_celltype'] <- S_data$celltype[i]}

identical(rownames(S_data),rownames(metadata[metadata$sub_celltype%in%S_data$celltype,]))#顺序和内容都一样

table(metadata$sub_celltype,metadata$celltype)

scRNA@meta.data <- metadata
table(scRNA$sub_celltype)
saveRDS(scRNA,file = 'data/scRNA_celltype.rds')

##check
rm(list = ls())
gc()
scRNA <- readRDS('data/scRNA_celltype.rds')
sc <- readRDS('data/scRNA.rds')
identical(sc$celltype,scRNA$celltype)
#TRUE

##plot
scRNA <- readRDS('data/scRNA_celltype.rds')
gene <- c('CD3E','CD8A','CD4','NCAM1','B3GAT1','TRDC','TRGC1','TRAC','TRBC1',#T/NK cell
          'FOXP3','IL2RA',#Treg
          'KLRG1','FGFBP2','CX3CR1','GZMB','PRF1','GNLY',#Temra         
          'CXCR5','CXCL13','PDCD1','BCL6','SOX4','PRDM1','ICOS','MME','TNFSF4','TNFRSF4','TNFRSF8',#TPH&TFH 
          'CD19','MS4A1','CD79A','PAX5','CD24','MME','CD27','IGHD','IGHE','IGHM','IGKC','IGLV',#B cell
          'CD38','SDC1','XBP1','PRDM1','CXCR4',#plasma
          'ITGAM','ITGAX','CD14','FCGR3A','CD68','CD163','FUT4','CD123',#mono/macro
          'CD1C','CLEC9A',#cDC
          'LILRA4','IRF7',#pDC
          'CR1','CR2','VCAM1',#FDC
          'FCGR3B','CEACAM8','S100A8','S100A9',#中性粒
          'PECAM1','VWF','CLDN5',#内皮         
          'MKI67'
          )
library(Seurat)
library(ggplot2)
Idents(scRNA) <- scRNA$sub_celltype
DotPlot(scRNA, features = unique(gene),
        assay='RNA'  )  + coord_flip()+
  scale_color_gradientn(colours = viridis::viridis(20))+
  theme(axis.text.x=element_text(angle = 45,vjust=0.5))
ggsave('result/celltype_dotplot.pdf',width = 10,height = 12)
table(scRNA$sub_celltype,scRNA$celltype)

####################################
##heatmap
####################################
avg_exp <- AverageExpression(scRNA, features = gene, group.by = "sub_celltype",assays = 'RNA',slot = 'data')
exp_matrix <- avg_exp$RNA

head(exp_matrix)
colnames(exp_matrix)
rownames(exp_matrix)

#group
sampledata <- data.frame(row.names = colnames(exp_matrix),
                         group = 1:ncol(exp_matrix))
rownames(sampledata) <- gsub('_','-',rownames(sampledata))#前面提取平均表达基因时系统自动将'_'替换为了'-',因此这里也要一样，防止识别错误
sampledata[rownames(sampledata) %in% c('gγδ-T','TFH1','TFH2','DNT','TPH','TEMRA-CD8-T','Treg','Ki67+-γδ-T','TFH3','DPT','TFH3'),1]='T_cell'
sampledata[rownames(sampledata) %in% c('naïve B','Plasma','Plasmablast','T-B'),1]='B_cell'
sampledata[rownames(sampledata) %in% c('pDC','Monocytes','Neutrophils','Macrophage(M2)','cDC'),1]='myeliod'
sampledata[rownames(sampledata) %in% c('FDC','endothelial-cells','T-FDC','undefined'),1]='stromal'
table(sampledata)

genedata <- data.frame(row.names = rownames(exp_matrix),
                         group = 1:nrow(exp_matrix))
genedata[rownames(genedata) %in% c('CD3E','CD8A','CD4','NCAM1','B3GAT1','TRDC','TRGC1','TRAC','TRBC1'),1]='T/NK_cell'
genedata[rownames(genedata) %in% c('FOXP3','IL2RA'),1]='Treg'
genedata[rownames(genedata) %in% c('KLRG1','FGFBP2','CX3CR1','GZMB','PRF1','GNLY'),1]='Temra'
genedata[rownames(genedata) %in% c('CXCR5','CXCL13','PDCD1','BCL6','SOX4','PRDM1','ICOS','MME','TNFSF4','TNFRSF4','TNFRSF8'),1]='TPH&TFH'
genedata[rownames(genedata) %in% c('CD19','MS4A1','CD79A','PAX5','CD24','MME','CD27','IGHD','IGHE','IGHM','IGKC','IGLV'),1]='B_cell'
genedata[rownames(genedata) %in% c('CD38','SDC1','XBP1','PRDM1','CXCR4'),1]='plasma'
genedata[rownames(genedata) %in% c('ITGAM','ITGAX','CD14','FCGR3A','CD68','CD163','FUT4','CD123'),1]='mono/macro'
genedata[rownames(genedata) %in% c('CD1C','CLEC9A'),1]='cDC'
genedata[rownames(genedata) %in% c('LILRA4','IRF7'),1]='pDC'
genedata[rownames(genedata) %in% c('CR1','CR2','VCAM1'),1]='FDC'
genedata[rownames(genedata) %in% c('FCGR3B','CEACAM8','S100A8','S100A9'),1]='Neutrophil'
genedata[rownames(genedata) %in% c('PECAM1','VWF','CLDN5'),1]='Endothelial_cells'
genedata[rownames(genedata) %in% c('MKI67'),1]='Proliferate_cells'
table(genedata)
identical(rownames(genedata),rownames(exp_matrix))
identical(rownames(sampledata),colnames(exp_matrix))
#plot
table(sampledata)
sampledata <- sampledata[order(sampledata$group),,drop=F]
exp_matrix <- exp_matrix[,match(rownames(sampledata),colnames(exp_matrix))]

genedata <- genedata[order(genedata$group),,drop=F]
exp_matrix <- exp_matrix[match(rownames(genedata),rownames(exp_matrix)),]
library(pheatmap)
pheatmap(
  exp_matrix,
  annotation_col = sampledata,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  fontsize_row = 8,
  fontsize_col = 6,
  gaps_row = c(11,13,16,19,26,30,32,37,38,47,53,62),
  gaps_col = c(4,9,13),
  angle_col = "45"
)


