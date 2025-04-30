

#1. 载入需要的工具包

rm(list = ls())#清除内存变量
library(flowCore)
library(Rcpp)
library(cytofkit)
library(igraph)
library(ggplot2)
library(ggthemes)
library(Rtsne)
library(dplyr)
library(cytofexplorer)

#2.数据预处理（读取->Downsample->合并->Transform）

#2.1读取FCS文件，Downsample，Merge等
#输入
projectdir<-"E:/work/AITL-CYTOF/clinical/Pipeline_1"
mergeMethod <- "ceil" #合并方法："ceil", "all", "fixed", "min"
fixedNum <- 10000       # 设置从每个文件抽取的细胞数目



#目录设置
wdir <-"21_PhenoGraph_tSNE"   #引号内更改为目标文件所在子目录
raw_fcs_dir="01_rawfcs"
wdir<-paste0(projectdir,"/",wdir)
raw_fcs_dir<-paste0(projectdir,"/",raw_fcs_dir)
metadata_dir<-paste0(projectdir,"/02_metadata")

setwd(wdir)

#读取文件，DownSample，Merge
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


#2.2 Transform 数据转换：仅针对部分通道，方法cytofAsinh
#input:
#重要：通道选择，仅仅首次需要input
#write.csv(colnames(combined_data_raw),paste0(metadata_dir,"/all_markers.csv"),row.names = FALSE)  #仅第一次运行，通道一旦选好可以将此行注释

#到工作目录中找到all_markers.csv,打开后，首行A列标明“markers”，B列“transform”
#B列中在要transform的标签行标记1，保存退出
#后面可以把write.csv.....这一行注释，就不会再运行此行


#数据转换
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
transform_id=which(all_markers$transform==1)
combined_data_transformed<-combined_data_raw
combined_data_transformed[, transform_id] <- apply(combined_data_transformed[, transform_id,drop = FALSE],
                                                   2,cytofAsinh)
#checkpoint
head(combined_data_transformed)




#3. 运行 PhenoGraph聚类(仅有一个参数K，默认30)

#输入
#PhenoGraph参数设置：
k=20  #计算Knn网络的“邻居”个数

#通道选择：到metadata目录中找到all_markers.csv,打开后，在右面选择一个空列输入“PhenoGraph”，
#并在要进行PhenoGraph聚类的标签行标记1，保存退出。


#模块主体部分
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
PhenoGraph_id=which(all_markers$PhenoGraph==1)
PhenoGraph_input_data=combined_data_transformed[,PhenoGraph_id]

#phenograph的elbow测试，由于要进行多轮聚类测试，耗时较长，建议在total event较少的情况下使用
#PG_elbow(PhenoGraph_input_data)

PhenoGraph_result <-as.numeric(membership(Rphenograph(data = PhenoGraph_input_data,k=k)))


#Checkpoint
hist(PhenoGraph_result,unique(PhenoGraph_result))





#4. 运行tSNE(BH-sne)降维(two options，二选一，本步骤需要时间较长)
#Input: tSNE参数设置：

max_iter=1000   #迭代次数
perplexity=30  #困惑度
seed=123      #随机数种子
theta=0.5      #权衡速度与准确度，越小越精确，越大速度越快
dims = 2       #降维输出维度（默认是2）

#通道选择：找到all_markers.csv,打开后，在新的一列首行输入“tSNE”，
#并在要进行tSNE降维的标签行标记1，保存退出。

#tsne分析:
if (exists('seed')) set.seed(seed)
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
tSNE_para_id=which(all_markers$tSNE==1)
tSNE_input_data=combined_data_transformed[,tSNE_para_id]

tsne_result <- Rtsne(tSNE_input_data,
                     initial_dims = ncol(tSNE_input_data),
                     pca = FALSE,
                     dims = dims,
                     check_duplicates = FALSE,
                     perplexity=perplexity,
                     max_iter=max_iter,
                     theta=theta)$Y
row.names(tsne_result)<-row.names(combined_data_raw)
colnames(tsne_result)<-c("tsne_1","tsne_2")

#checkpoint: 以下语句可以绘制出降维分析的草图
plot(tsne_result)





#5降维和聚类结果可视化
#5.1：降维和聚类结果做图（Shinny App）

#input:
##降维结果汇总
dimReducedRes <- list(tsne_result)  #降维分析结果（可以将多次结果汇总）
names(dimReducedRes)<-c("tsne")     #降维分析名称（可以将多次结果汇总，名称与前面结果顺序一致）
str(dimReducedRes)

clusterRes<-list(PhenoGraph_result) #聚类分析结果（可以将多次结果汇总）
names(clusterRes)<-c("Phenograph")  #聚类分析名称（可以将多次结果汇总，名称与前面结果顺序一致）
str(clusterRes)

source("./backup/export_cytofkit_RData.R")

#生成Cytofkit_ShinnyApp 识别的RData文件
export_cytofkit_RData(expressionData = combined_data_raw,
                      dimReducedRes = dimReducedRes,
                      clusterRes = clusterRes,
                      projectName = "cytofkit",
                      rawFCSdir = raw_fcs_dir,
                      resultDir = raw_fcs_dir,
                      dimRedMarkers = colnames(tSNE_input_data),
                      sampleNames = unique(File_ID))

#如看结果，请新打开一个Rstudio窗口，打开Open_ShinnyApp.R 运行里面的语句。点击YES，在浏览器操作。


#追加MetaData数据

#确定分组
#文件分组，仅仅首次需要input
#write.csv(unique(File_ID),paste0(metadata_dir,"/all_samples.csv"),row.names = FALSE)

#'到metadata目录中找到samplename.csv,打开后，
#'A列名设为File_ID，根据样本分组情况设置后面几列；
#'B列明设为Short_name，每个样本的简称，可以方便的出现在输出的图片中；
#'C列以后根据实验设计，输入实际的分组情况例如：Timepiont, Tissue_Type,Patient_ID ect.
#'第一次运行，通道一旦选好可以将此行注释

groups<-read.csv(paste0(metadata_dir,"/all_samples.csv"),header = TRUE,stringsAsFactors=FALSE)
groups$File_ID<-unique(File_ID)

#checkpoint
groups

#5.2 降维和聚类结果做图（R语言脚本直接生成）
#13. 整合聚类和降维数据
source("./backup/FlowSOM_tSNE_Output.R")
combined_data_plot <- cbind(combined_data_transformed,
                                tsne_result,
                                PhenoGraph = PhenoGraph_result)

#checkpoint
head(combined_data_plot)
#source("./backup/FlowSOM_tSNE Analysis backup (Pipeline2-2).R")
#生成tsne-Phenograph系列图片
draw_tsne_figs(combined_data_plot=combined_data_plot,
               groups=groups,
               cluster_color=dif_seq_rainbow,
               cluster_name="PhenoGraph",
               major_cond="group1",#这里要改成分组名
               reduction_dm1="tsne_1",
               reduction_dm2="tsne_2")

#生成Marker-Cluster 系列Density图片, 在。csv中加一列“density_plot”，想要对比的marker标1
all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
draw_density_plots(combined_data_plot,
                   groups=groups,
                   all_markers=all_markers,
                   cluster_color=dif_seq_rainbow,  #选择颜色 其他选择：rainbow，dif_seq_rainbow,brewer_color_sets
                   cluster_name="PhenoGraph",      #cluster 通道的名称
                   cluster_id=c(20, 31)                # 选择要显示的cluster
                   )

#生成所有marker的tsne热图
#'打开all_markers.csv，在后面的一个空列首行标记”tsne_heatmap“,
#'然后在要生成热图的所在行标记1，保存退出

all_markers=read.csv(paste0(metadata_dir,"/all_markers.csv"))
all_markers$markers<-colnames(combined_data_raw)[1:nrow(all_markers)]
heatmap_tsne_id<-which(all_markers$tsne_heatmap==1)
heatmap_tsne_markers=colnames(combined_data_transformed[,heatmap_tsne_id])
draw_tsne_heatmaps(combined_data_plot=combined_data_plot,
                   heatmap_tsne_markers=heatmap_tsne_markers,
                   single_file=T)





#6 将结果导出成FCS文件

## Input： 整合待导出数据
combined_data_output <- data.frame( tsne_result,
                                    PhenoGraph = PhenoGraph_result)
#head(combined_data_output)

#根据File_ID将合并数据还原成单个样本数据，导出FCS文件
row.names(combined_data_output)<-row.names(combined_data_raw)
cytof_addToFCS_modified(
  data=combined_data_output,
  rawFCSdir=raw_fcs_dir,
  analyzedFCSdir=paste(wdir,"/output",sep = ""),
  newHeader = "pg_tsne_"
)

