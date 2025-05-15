rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(ggplot2)

data.path <- 'raw-data'
data.list <- list.files(data.path)
data.list

scRNA_list <- list()
for (i in 1:length(data.list)){
  print(data.list[i])
  path <- paste0(data.path,'/',data.list[i])
  data <- Read10X(path)
  seurat_data <- CreateSeuratObject(counts = data,min.features = 200,min.cells = 3,project = data.list[i])
  scRNA_list[[i]] <- seurat_data
}
####################################
#QC
####################################
dir.create('QC')
for(i in 1:length(scRNA_list)){
  sc <- scRNA_list[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)]
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)
  scRNA_list[[i]] <- sc
  rm(sc)
}

violin_before <- list()
for(i in 1:length(scRNA_list)){
  violin_before[[i]] <- VlnPlot(scRNA_list[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}
violin_before

for(i in 1:length(violin_before)){
  print(violin_before[[i]])
  ggsave(filename = paste0('QC/',levels(violin_before[[i]][[1]][["data"]][["ident"]]),'QC_before.pdf'),width = 10,height = 10)
}


scRNA_list <- lapply(X = scRNA_list, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & 
                mt_percent < 20 & 
                HB_percent < 5 & 
                nCount_RNA < quantile(nCount_RNA,0.97))})

scRNA_list <- merge(x=scRNA_list[[1]],y=scRNA_list[-1])
violin_after <- VlnPlot(scRNA_list,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                        pt.size = 0.01,
                        ncol = 4)
violin_after
ggsave('QC/QC_after.pdf',width = 15,height = 10)

dim(scRNA_list)
####################################
#Data normalization and feature selection
####################################
scRNA_list <- NormalizeData(scRNA_list) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = T)
a=DimPlot(scRNA_list,reduction = "pca")
a

scRNA_list <- JoinLayers(scRNA_list)

dir.create('harmony')
ggsave('harmony/harmony_before.pdf',height = 5,width = 8)
dir.create('data')
saveRDS(scRNA_list,file = 'data/scRNA.rds')
####################################
#harmony
####################################
rm(list = ls())
gc()
scRNA <- readRDS('data/scRNA.rds')
library(harmony)
#BiocManager::install('harmony')
table(scRNA$orig.ident)
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident")
scRNA@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA,reduction = "harmony",group.by = "orig.ident")
b
ggsave(filename = 'harmony/harmony_after.pdf',height = 5,width = 8)

ElbowPlot(scRNA, ndims=50, reduction="harmony")
pct <- scRNA[["pca"]]@stdev / sum(scRNA[["pca"]]@stdev) * 100 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use
####################################
#FindNeighbors
####################################
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = seq(0.4,1,by=0.1))
scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:20)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)

saveRDS(scRNA,file = 'data/scRNA.rds')
