
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
