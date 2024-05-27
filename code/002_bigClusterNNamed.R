suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(ggplot2)
    library(SCP)
    library(cowplot)
    library(RColorBrewer)
    library(patchwork)
    library(ggpubr)
    library(qs)
})

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  读入聚类结果，使用singleR和经典marker基因进行注释
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 读入数据聚类的结果***********************************************************************************************
setwd("/root/wangje/Project/OveryArtical")
scRNA <- qread("/root/wangje/Project/OveryArtical/04_大群数据_scp.qs")
# 指定reduction
reduction <- "HarmonyUMAP2D"
# 绘制该reduction下不同分辨率的结果
Seurat::DimPlot(scRNA, reduction = reduction, group.by = paste0("Harmony_SNN_res.", seq(0.1, 1, 0.1)), ncol = 5, raster = T, label = T, repel = T) %>%
    ggsave(filename = file.path(paste0(reduction, "大群不同分辨率聚类结果.png")), height = 8, width = 25, bg = "white")

# 添加信息**********************************************************************************************************
Seurat::DimPlot(scRNA, reduction = reduction, group.by = c("sample", "group"), ncol = 2, raster = T, label = F, repel = F) %>%
    ggsave(filename = file.path(paste0(reduction, "大群batch结果.png")), height = 4, width = 10, bg = "white")
# 添加年龄
# unique(scRNA$sample)
#  [1] "HRS421447" "HRS421448" "HRS421449" "HRS421450" "HRS421451" "HRS421452" "HRS421453" "young1"    "young3"    "young4"
# [11] "middle2"   "middle3"   "middle4"   "old1"      "old2"      "old3"      "na-1010"   "na-426"    "na1111"    "na1122"
# [21] "na1117"    "na1128"    "NA_728"    "na-412"    "na1129"
library(magrittr)
scRNA %<>% dplyr::mutate(
    Age = dplyr::case_when(
        sample %in% c()
    )
)

# 使用singleR进行粗注释***********************************************************************************************
suppressMessages({
    library(Seurat)
    library(dplyr)
    library(patchwork)
    library(ggplot2)
    library(SingleR)
})

#### singleR ####
Run_singleR <- function(data, cluster, resolution) {
    require("SingleR")
    require("dplyr")
    stopifnot(file.exists("/root/wangje/singleR.RData"))
    load("/root/wangje/singleR.RData")
    ref_list <- list(encode, hema, hpca, immune, monaImmune)
    labels_list <- list(
        encode$label.main,
        hema$label.main,
        hpca$label.main,
        immune$label.main,
        monaImmune$label.main
    )
    sce.data <- GetAssayData(data, solt = "data")
    sce.singleR <- SingleR::SingleR(
        test = sce.data,
        ref = ref_list,
        labels = labels_list,
        clusters = cluster,
        BPPARAM = BiocParallel::MulticoreParam(6)
    )
    # 提取 celltype 数据
    celltype <- data.frame(
        ClusterID = rownames(sce.singleR),
        singleR = sce.singleR$labels,
        stringsAsFactors = FALSE
    )
    # 将注释信息添加到seurat对象中
    if (is.null(resolution)) stop("missing resolution files.")
    tmp <- data@meta.data %>% dplyr::select(resolution)
    colnames(tmp) <- "ClusterID"
    tmp1 <- left_join(tmp, celltype, by = "ClusterID")
    rownames(tmp1) <- rownames(tmp)
    sce_new <- AddMetaData(data, tmp1)
    return(sce_new)
}

scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$Harmony_SNN_res.1,
    resolution = "Harmony_SNN_res.1"
)
scRNA$Harmony_singleR <- scRNA$singleR # 列名更改名称

# 查看注释后的结果*************************************************************************************************
plotBigDotPlot <- function(inputfile, group.by, marker) {
    p1 <- DotPlot(inputfile, assay = "RNA", features = marker, scale = TRUE, group.by = group.by) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red") +
        # 添加四条黑色边界线
        annotate(geom = "segment", x = -Inf, xend = Inf, y = Inf, yend = Inf, color = "black", size = 1) +
        annotate(geom = "segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, color = "black", size = 1) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
            strip.text.x = element_text(size = 14, angle = 0),
            axis.text.y = element_text(size = 15),
            panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7, linetype = "solid"),
            legend.position = "bottom",
            panel.spacing.x = unit(0.1, "cm"),
            panel.spacing.y = unit(0.1, "cm"),
            panel.grid = element_line(color = "grey", linetype = "dashed", linewidth = unit(0.1, "cm"))
        ) +
        labs(x = "Gene Marker", y = " ")
    return(p1)
}

# 绘制marker基因表达图*******************************************************************************************
mark_list <- list(
    Immune = c("PTPRC"),
    NK = c("KLRC1", "FCGR3A", "NCAM1", "KLRD1", "GNLY"),
    "T cells" = c("CD3D", "CD3G", "CD2", "CD8A", "CD8B", "CD4"),
    Fibroblasts = c("COL1A1", "DCN", "LUM", "PDGFRA"),
    Myeloids = c("LYZ", "CD68", "TYROBP"),
    Epithelial = c("CD24", "KRT19", "EPCAM", "KRT18"),
    Bcells = c("CD79A", "CD19", "MS4A1"),
    Endothelial = c("CLDN5", "FLT1", "RAMP2", "CDH5"),
    "Smooth muscle" = c("GJA4", "PLAC9", "CRYAB", "ACTA2", "PLN", "ADIRF", "MYH11", "CNN1", "MYL9"),
    Neutrophil = c("FCGR3B", "S100P", "CMTM2"),
    oocytes = c("TUBB8", "ZP3", "FIGLA"),
    Granulosal = c("GSTA1", "AMH", "HSD17B1", "DSP"),
    DC = c("LILRA4", "CXCR3", "IRF7"),
    Mast = c("CPA3", "TPSAB1", "TPSB2")
)
p1 <- Seurat::DimPlot(scRNA, group.by = "singleR", label = T, repel = T, raster = F, reduction = reduction) + labs(xlab = "umap1", ylab = "umap2")
p2 <- plotBigDotPlot(scRNA, group.by = "singleR", marker = mark_list)

p3 <- Seurat::DimPlot(scRNA, group.by = "Harmony_SNN_res.1", label = T, repel = T, raster = F, reduction = reduction) + labs(xlab = "umap1", ylab = "umap2")
p4 <- plotBigDotPlot(scRNA, group.by = "Harmony_SNN_res.1", marker = mark_list)
# 保存图片
ggsave(
    filename = file.path(paste0(reduction, "_singleR注释结果.png")), height = 17, width = 34, limitsize = F,
    plot = cowplot::plot_grid(plotlist = list(p1, p2, p3, p4), align = "hv", ncol = 2, nrow = 2, rel_widths = c(1.3, 4)), bg = "white"
)

# 添加celltype注释信息*******************************************************************************************
library(magrittr)
scRNA@meta.data %<>% dplyr::mutate(
    celltype = dplyr::case_when(
        
    )
)