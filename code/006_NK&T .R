
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

# 1 提取NK和T細胞亚群进行聚类分析**************************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA,subset=celltype == "NK & T cells")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,meta.data = scRNA@meta.data[,c('sample')])
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'Harmony',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './07_T细胞数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'BBKNN',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './07_T细胞数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'scVI',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './07_T细胞数据_scp.qs')

# 2 使用singleR粗注释*************************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$BBKNN_res.1.5,
    resolution = "BBKNN_res.1.5"
)

# 3 使用NK 和T细胞的marker进行注释**************************************************************************************************
markers <- function(){
    marks <- list(
        TCells = c('CD3G','CD3D','CD2','CD4','CD8A','CD8B'),
        NK = c('FGFBP2','FCG3RA','CX3CR1','CD56','CXCR3','IFNG','KLRC1','IL7R','KIT','KLRD1','GZMB'))}


Idents(scRNA) <- scRNA$BBKNN_res.1.5
p1 <- plotBigDotPlot(scRNA,group.by = 'BBKNN_res.1.5',marker = markers())