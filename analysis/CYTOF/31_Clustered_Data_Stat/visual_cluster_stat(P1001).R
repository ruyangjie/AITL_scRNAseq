
#1. 载入需要的工具包

rm(list = ls())#清除内存变量
gc()
library(cytofkit)
library(dplyr)
library(gplots)
library(ggrepel)
library(Rmisc)
library(colorRamps)
library(RColorBrewer)
library(car)
library(ggpubr)
library(cytofexplorer)
library(reshape2)


#2.数据预处理（读取->Downsample->合并->Transform）

#2.1读取->Downsample->合并
projectdir<-"E:/work/AITL-CYTOF/clinical/Pipeline_1"
mergeMethod <- "all" #合并方法："ceil", "all", "fixed", "min"
fixedNum <- 10000         # 设置从每个文件抽取的细胞数目



#目录设置
wdir <-"31_Clustered_Data_Stat"   #引号内更改为目标文件所在子目录
raw_fcs_dir="21_PhenoGraph_tSNE/output"

#读取文件，DownSample，Merge
wdir<-paste0(projectdir,"/",wdir)
raw_fcs_dir<-paste0(projectdir,"/",raw_fcs_dir)
metadata_dir<-paste0(projectdir,"/02_metadata")
setwd(wdir)


file_name <- list.files(raw_fcs_dir,pattern='.fcs$', full=TRUE)
combined_data_raw <- cytof_exprsMerge(fcsFiles = file_name,
                                      transformMethod = "none",
                                      mergeMethod =mergeMethod,
                                      fixedNum = fixedNum)
#简化列名
paraname<-colnames(combined_data_raw)
paraname<-sub(".*<","",paraname)
paraname<-sub(">.*","",paraname)
paraname<-sub("-","_",paraname)
#paraname<-sub("+","_",paraname)
colnames(combined_data_raw)<-paraname

#增加File_ID
File_ID<-sub("_[0-9]*$","",row.names(combined_data_raw))
combined_data_raw<-data.frame(combined_data_raw,File_ID)


#Checkpoint
head(combined_data_raw)


#2.2 marker选项
#打开all_markers.csv，增加expr_para、heatmap两列分别用来指定用来分析差异表达和出现在Heatmap中的marker

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),header = TRUE)
all_markers$markers<-colnames(combined_data_raw)[1:nrow(all_markers)]

#checkpoint
all_markers


#2.3 确定分组
groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors=FALSE)
groups$File_ID<-unique(File_ID)

#checkpoint
groups



#3针对各个cluster进行统计分析
# 
source("./backup/Cluster_Expr.R")
source("./backup/Cluster_Abundance.R")

#3.1#数据初步整理，全局参数设定
cluster_stat<- stat_by_cluster( combined_data_raw,
                                all_markers,
                                cluster_name="PhenoGraph",
                                summerise_method="median",
                                groups=groups,
                                major_cond="group3",
                                group_seq=c("0","1"), #设置major_cod 在各个统计图中的顺序
                                stat.paired=F,
                                stat.method="t-test"
                                )

#3.2亚群丰度
#绘制boxplot
draw_abundance_boxplot(cluster_stat,
                       #cluster_id=c(1,2,3,5),
                       boxplot_ctrl="0"
                       #comparisons = list(c("PBMC","Biopsy"))
                       )


#绘制亚群丰度火山图
draw_abundance_volcano(cluster_stat,
                       cond1="0",
                       cond2="1")


#3.3 Marker差异表达分析
##需要all_markers.csv中含有expr_para一列
#生成每个cluster的expr_para的heatmap,boxplot
cluster_expr_report(cluster_stat,
                    cluster_id=c(1,10,22),  
                    heatmap_ctrl="NO"
                    #comparisons = list(c("PBMC","Biopsy"))
                    )


#绘制差异表达火山图
draw_expr_volcano(cluster_stat,
                  cond1="NO",
                  cond2="YES")



#3.4 绘制cluster的heatmap
##需要all_markers.csv中含有heatmap一列
source("./backup/Heatmap_Output.R")

cluster_expr_matrix<-cluster_stat$cluster_expr_matrix

#'到metadata目录中找到samplename.csv,打开后，
#'在后面添加一列，列名为heatmap，所有需要在热图上显示的marker上标记1
#'多余的列可以保留，只要列名不重复，就不会影响程序运行。
#'保存后退出。

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"),header = TRUE)
heatmap_ID<-as.character(dplyr::filter(all_markers,heatmap==1)$markers)
draw_expr_heatmap(xdata=cluster_expr_matrix[,heatmap_ID],
                  Rowv=T,
                  Colv=T,
                  trans_method="CytofAsinh",
                  color_style=0,
                  colorkeys = rev(brewer.pal(n = 7, name ="RdYlBu"))
                  )


