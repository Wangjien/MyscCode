

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

# 1 提取卵细胞和颗粒细胞亚群进行聚类分析********************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA, subset = celltype == "Oocyte|granule cells")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts)
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data), "_[A|T|G|C].*", n = 2)[, 1]
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./08_内皮胞数据_scp.qs")
# 2 使用singleR粗注释*************************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$Harmony_SNN_res.1,
    resolution = "Harmony_SNN_res.1"
)