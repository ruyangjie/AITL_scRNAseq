library(Seurat)
library(ggrepel)
library(dplyr)
rm(list = ls())

scRNA <- readRDS('data/scRNA_celltype.rds')
Idents(scRNA) <- 'sub_celltype'

#Differential analysis
table(scRNA$sub_celltype)
markers = FindMarkers(scRNA,
                        logfc.threshold = 0,
                        min.pct = 0.1,
                        only.pos = FALSE,
                        ident.1 = "TPH",ident.2 = "TFH2")%>%
  mutate(gene = rownames(.))

getwd()
dir.create('result')
write.csv(markers,file = "result/TPHvsTFH(control).csv")

####################################
#plot
####################################
rm(list = ls())
data <- read.csv("result/TPHvsTFH(control).csv",row.names = 1)
data %>%
  mutate(change = case_when(
    avg_log2FC >= 0 & p_val_adj <=0.05 ~ "UP",
    avg_log2FC <= 0 & p_val_adj <=0.05 ~ "DOWN",
    TRUE ~ "NONE"
  )) -> data
table(data$change)
#markergene
gene <- c('IL21','CD70','PTGES2','PTGES3','STAT5A','NR4A2','SOCS1','TGFB1','MIF','BCL6','EMP3')
data$marker_gene <- NA
data$marker_gene[match(gene,rownames(data))] <- data$gene[match(gene,rownames(data))]


ggplot(data,aes(x = avg_log2FC,y= -log10(p_val_adj)))+
  geom_point(aes(color=change),size= 1)+
  ggtitle(label = "TPH_vs_TFH2")+
  theme_bw()+
  theme_classic()+
  theme(
    plot.title = element_text(face="bold",hjust = 0.5),
    plot.background = element_rect(fill = "transparent",color = NA)
  )+
  scale_x_continuous(limits = c(-7, 7))+
  scale_y_continuous(limits = c(0, 300))+
  geom_vline(xintercept = c(0),lty = 2, col = "black", lwd = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5)+
  scale_colour_manual(values = c('#00B2BF','#898989','#C82E31'))
ggsave("TPH_vs_TFH2_火山图_p_adj.pdf",width = 8,height = 7)

####################################
#p_value
####################################
rm(list = ls())
gc()
data <- read.csv("result/TPHvsTFH(control).csv",row.names = 1)

data %>%
  mutate(change = case_when(
    avg_log2FC >= 0 & p_val <=0.05 ~ "UP",
    avg_log2FC <= 0 & p_val <=0.05 ~ "DOWN",
    TRUE ~ "NONE"
  )) -> data
table(data$change)
#markergene
gene <- c('IL21','CD70','PTGES2','PTGES3','STAT5A','NR4A2','SOCS1','TGFB1','MIF','BCL6','EMP3')

data$marker_gene <- ''
data$marker_gene[match(gene,rownames(data),nomatch = NA)] <- data$gene[match(gene,rownames(data),nomatch = NA)]


ggplot(data,aes(x = avg_log2FC,y= -log10(p_val)))+
  geom_point(aes(color=change),size= 1)+
  ggtitle(label = "TPH_vs_TFH2")+
  geom_text_repel(aes(label=marker_gene),size=3,color="#0023F5",max.overlaps =10000,nudge_x = 0.5,nudge_y = 1)+#,
  theme_bw()+
  theme_classic()+
  theme(
    plot.title = element_text(face="bold",hjust = 0.5),
    plot.background = element_rect(fill = "transparent",color = NA)
  )+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(0, 100))+
  geom_vline(xintercept = c(0),lty = 2, col = "black", lwd = 0.5)+
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.5)+
  scale_colour_manual(values = c('#00B2BF','#898989','#C82E31'))
ggsave("TPH_vs_TFH2_火山图_p_val.pdf",width = 12,height = 10)

####################################
##tph_vs_other_cluster
####################################
rm(list = ls())
scRNA <- readRDS('data/scRNA_celltype.rds')
Idents(scRNA) <- 'sub_celltype'
table(scRNA$sub_celltype)
markers = FindAllMarkers(scRNA,
                         logfc.threshold = 0,
                         min.pct = 0.01,
                         only.pos = FALSE)%>%
  mutate(gene = rownames(.))
write.csv(markers,'find_allmarker.csv')
