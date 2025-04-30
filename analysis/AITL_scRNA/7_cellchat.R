rm(list = ls())
gc()
#devtools::install_github("jinworks/CellChat")
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(magrittr)
library(CellChat)

#devtools::install_github('immunogenomics/presto')
#prepare
data <- readRDS('T_B_cellchat/cellchat_use_data.rds')
Idents(data) <- 'sub_celltype'
table(data$sub_celltype)
data.object <- data
data.input <- GetAssayData(data.object, assay = "RNA", layer = "data")
data.meta <- data.object@meta.data[,c('celltype','sub_celltype')]

#cellchat
data.cellchat <- createCellChat(object = data.input)
data.cellchat <- addMeta(data.cellchat, meta = data.meta)
data.cellchat <- setIdent(data.cellchat, ident.use = "sub_celltype")
table(data.cellchat@idents)

CellChatDB <- CellChatDB.human
data.cellchat@DB <- CellChatDB

options(future.globals.maxSize = 3e9)
data.cellchat <- subsetData(data.cellchat, features = NULL)
options(future.seed = T)

data.cellchat <- identifyOverExpressedGenes(data.cellchat)
data.cellchat <- identifyOverExpressedInteractions(data.cellchat)
data.cellchat <- smoothData(data.cellchat, adj = PPI.human)


data.cellchat <- computeCommunProb(data.cellchat,raw.use = T)



data.cellchat <- filterCommunication(data.cellchat,min.cells = 10)

data.cellchat <- computeCommunProbPathway(data.cellchat)

data.cellchat <- aggregateNet(data.cellchat)


data.cellchat <- netAnalysis_computeCentrality(data.cellchat,slot.name = "netP")


group.net <- subsetCommunication(data.cellchat)

write.csv(group.net,file = "T_B_cellchat/group_net_inter_raw.use.csv")
saveRDS(data.cellchat,file = 'T_B_cellchat/cellchatdata.rds')












#plot

rm(list = ls())
gc()
data.cellchat <- readRDS('T_B_cellchat/cellchatdata.rds')
groupSize <- as.numeric(table(data.cellchat@idents))
par(mfrow = c(1,2))
netVisual_circle(data.cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 #arrow.size = 0.1,
                 title.name = "Number of interactions")


netVisual_circle(data.cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

mat <- data.cellchat@net$weight

par(mfrow = c(4,5), mar = c(1,1,1,1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize,
                   arrow.width = 0.2,arrow.size = 0.05,
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}


#plot2
data.cellchat@netP$pathways
pathways.show <- 'SELPLG'
summary(data.cellchat@net)
par(mfrow=c(1,1))

netVisual_aggregate(data.cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")

pathways.show <- 'CD70'
par(mfrow=c(1,1))
netVisual_aggregate(data.cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")

pathways.show <- 'Prostaglandin'
par(mfrow=c(1,1))
netVisual_aggregate(data.cellchat, 
                    signaling = pathways.show, 
                    layout = "circle")

#SELPLG
pathways.show <- 'SELPLG'#IL21,SELPLG,CD70
LR <- extractEnrichedLR(data.cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR
LR.show <- LR[LR$interaction_name=='SELPLG_SELL',]

netVisual_individual(data.cellchat, signaling = pathways.show,  pairLR.use = LR.show,layout = "circle")
netVisual_heatmap(data.cellchat, signaling = pathways.show, color.heatmap = "Reds")
