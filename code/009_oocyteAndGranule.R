
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
qsave(scRNA, file = "./10_卵细胞和颗粒细胞数据_scp.qs")
# 2 使用singleR粗注释*************************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$Harmony_SNN_res.1,
    resolution = "Harmony_SNN_res.1"
)
DimPlot(scRNA, group.by = "singleR", reduction = "HarmonyUMAP2D", label = T, repel = T) %>%
    ggsave(filename = "/10_卵细胞和颗粒细胞singleR粗注释结果.png", height = 4, width = 5, bg = "white")
qsave(scRNA[, scRNA$singleR == "Fibroblasts"], file = "10_卵细胞和颗粒筛选出成纤维细胞.qs")
qsave(scRNA[, scRNA$singleR == "CD8+ T-cells"], file = "10_卵细胞和颗粒筛选T细胞.qs")
qsave(scRNA[, scRNA$singleR == "Epithelial cells"], file = "10_卵细胞和颗粒筛选上皮细胞.qs")
qsave(scRNA[, scRNA$singleR %in% c("Endothelial cells")], file = "10_卵细胞和颗粒筛选内皮细胞_scp.qs")
qsave(scRNA[,scRNA$singleR %in% c('Gametocytes','iPS_cells','Mesangial cells','MSC,Neurons')],file = '10_卵细胞和颗粒细胞数据_scp.qs')