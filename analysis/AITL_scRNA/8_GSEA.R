library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
rm(list = ls())
gc()

gene <- read.csv("result/TFH2vsTPH(control).csv")
ids.data <- bitr(gene$X,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = "org.Hs.eg.db")
data=merge(gene,ids.data,by.x="X",by.y="SYMBOL")

data <- data[order(data$avg_log2FC,decreasing = T),]

cluster_GSEA_data <- as.numeric(data$avg_log2FC)
names(cluster_GSEA_data) <- data$ENTREZID
head(cluster_GSEA_data)

cluster_GSEA <- gseKEGG(cluster_GSEA_data,
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        nPerm = 10000,
                        pvalueCutoff = 0.05,
                        minGSSize = 10,
                        maxGSSize = 500)

GSEA_result <- cluster_GSEA@result
write.csv(GSEA_result,file = 'GSEA/TFH2vsTPH(control)_GSEA_result.csv')
saveRDS(cluster_GSEA,file = 'GSEA/TFH2vsTPH(control)_cluster_GSEA.rds')

#plot
cluster_GSEA <- readRDS('GSEA/TFH2vsTPH(control)_cluster_GSEA.rds')
pathway <- c('hsa04110','hsa03030','hsa03410','hsa03430')
color <- c("#f7ca64","#43a5bf","#86c697","#ef998a")#颜色选择


gseaplot2(cluster_GSEA,
          color = color,
          pvalue_table = F,
          geneSetID = pathway,
          base_size = 14)




gse_plot <- list()
for(i in pathway){

  plot <- gseaplot2(cluster_GSEA,
                    color = color,
                    pvalue_table = F,
                    geneSetID = i,
                    base_size = 14)
  gse_plot[[i]] <- plot
}
pdf(paste0('GSEA/result/GSEA_.pdf'),width = 5,height = 5,)
gse_plot
dev.off()
