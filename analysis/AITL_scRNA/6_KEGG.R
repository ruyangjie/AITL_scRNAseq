rm(list = ls())
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
data <- read.csv("result/TFH2vsTPH(control).csv")
data %>%
  mutate(change = case_when(
    avg_log2FC >= 1 & p_val <=0.05 ~ "UP",
    avg_log2FC <= -1 & p_val <=0.05 ~ "DOWN",
    TRUE ~ "NONE"
  )) -> data
table(data$change)

data <- data[data$change!="NONE",]
ids=bitr(data$X,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = "org.Hs.eg.db")
data.filter=merge(data,ids,by.x="X",by.y="SYMBOL")
data.id <- unique(data.filter$ENTREZID)
table(duplicated(data.id))
cluster <- data.id

#KEGG
dir.create('KEGG')
cluster_KEGG <- enrichKEGG(gene = cluster,organism = "hsa",pvalueCutoff = 0.05)
head(cluster_KEGG)
dotplot(cluster_KEGG,showCategory=10, title="TFH2 vs TPH KEGG")
ggsave(filename = 'result/TFH2_vs_TPH_KEGG.pdf',width = 6,height = 5)

write.csv(cluster_KEGG@result,file = "KEGG/TFH2vsTPH(control)_KEGG.CSV")
rm(list = ls())


cluster_KEGG <- read.csv("KEGG/TFH2vsTPH(control)_KEGG.CSV",row.names = 1)
data <- cluster_KEGG
cluster_KEGG <- cluster_KEGG[order(cluster_KEGG$Count,decreasing = T),]


data <- cluster_KEGG[1:10,]

ggplot(data,aes(x = Description,y = Count,width=0.6,fill = p.adjust))+
  geom_col()+
  coord_flip()+
  ggtitle('TFH2_vs_TPH_KEGG')+
  theme_bw()+
  scale_x_discrete(name = "KEGG PATHWAY")+
  scale_fill_gradient(low = "#D66460", high = "#3478AF")+ 
  theme(plot.title = element_text(hjust = 0.5,face = 'bold'),
        axis.title.x=element_text(hjust = 0.5,size = 15,face = 'bold'),
        axis.title.y=element_text(hjust = 0.5,size = 10,face = 'bold'),
        axis.text.y = element_text(size = 10),
        panel.border = element_rect( size = 2, fill = NA),
        panel.grid = element_blank())
ggsave(filename = 'KEGG/TFH2_VS_TPH_kegg.pdf',height = 5,width = 8) 

