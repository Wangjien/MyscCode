
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
setwd("/root/wangje/Project/OveryArtical")
scRNA <- qread("04_大群数据_scp.qs")

# 1 提取内皮細胞亚群进行聚类分析**************************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA, subset = celltype == "Endothelial cells")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = scRNA@meta.data[, c("sample", "group")])
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data), "_[A|T|G|C].*", n = 2)[, 1]
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮胞数据_scp.qs")
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "BBKNN", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮细胞数据_scp.qs")
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "scVI", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮细胞数据_scp.qs")
# 2 使用singleR粗注释*************************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$Harmony_SNN_res.1,
    resolution = "Harmony_SNN_res.1"
)
# Adipocytes Endothelial cells Endothelial_cells       Fibroblasts
#     1159             17881              1122              1910
DimPlot(scRNA, group.by = "singleR", reduction = "HarmonyUMAP2D", label = T, repel = T) %>%
    ggsave(filename = "08_内皮细胞singleR粗注释结果.png", height = 4, width = 5, bg = "white")
qsave(scRNA[,scRNA$singleR=="Fibroblasts"],file = "08_内皮细胞数据筛选出成纤维细胞.qs")
qsave(scRNA[,scRNA$singleR=="Adipocytes"],file = "08_内皮细胞数据筛选出脂肪细胞.qs")
qsave(scRNA[,scRNA$singleR%in%c('Endothelial cells','Endothelial_cells')],file = '08_内皮细胞数据_scp.qs')
# 3 去除其他细胞亚群后再聚类********************************************************************************************************
scRNA <- qread('08_内皮细胞数据_scp.qs')
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = scRNA@meta.data[, c("sample", "group")])
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data), "_[A|T|G|C].*", n = 2)[, 1]
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮胞数据_scp.qs")
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "BBKNN", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮细胞数据_scp.qs")
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "scVI", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮细胞数据_scp.qs")