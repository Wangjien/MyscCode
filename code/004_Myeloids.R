
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(ggplot2)
    library(SCP)
    library(tidydr)
    library(cowplot)
    library(RColorBrewer)
    library(patchwork)
    library(ggpubr)
    library(qs)
})
# 读入数据
setwd('/root/wangje/Project/OveryArtical')
scRNA <- qread('04_大群数据_scp.qs')

# 1 提取髓系亚群进行聚类分析**************************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA,subset=celltype == "Myeloid cells")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,meta.data = scRNA@meta.data[,c('sample','group')])
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'Harmony',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'BBKNN',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'scVI',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')

# 2 使用singleR粗注释*********************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$BBKNN_res.1,
    resolution = "BBKNN_res.1"
)
table(scRNA$singleR)
# Adipocytes        Basophils               DC  Dendritic cells Epithelial cells      Macrophages        Monocytes      Neutrophils
#     1200              418               49              170               23             4090             4978             2341
# 分布保存数据
qsave(scRNA[,scRNA$singleR == "Adipocytes"],file = "06_髓系细胞中提取脂肪细胞.qs")
qsave(scRNA[,scRNA$singleR == "Epithelial cells"],file = "06_髓系细胞中提取上皮细胞.qs")
qsave(scRNA[,!scRNA$singleR %in%c('Adipocytes','Epithelial cells')],file = "06_髓系细胞.qs")
# 使用singleR进行粗注释，可以发现在提取的髓系亚群中是有其他细胞存在的，可以使用FeaturePlot查看注释的准确性
# 去除其他细胞后再次进行聚类
