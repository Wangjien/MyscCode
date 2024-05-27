
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

# 1 提取髓系亚群进行聚类分析*******************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA,subset=celltype %in% c( "Stromal & theca cells","SMC"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,meta.data = scRNA@meta.data[,c('sample','group')])
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data),'_[A|T|G|C].*',n=2)[,1]
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'Harmony',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_基质细胞数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'BBKNN',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06基质细胞数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'scVI',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_基质细胞数据_scp.qs')

