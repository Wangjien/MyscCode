suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(magrittr)
  library(ggpubr)
  library(stringr)
  library(qs)
}))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  合并两个文章的数据与自测数据，并保存raw Seurat对象
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#**********************************************************************************************
# Single-Cell Atlas of Human Ovaries Reveals The Role Of The Pyroptotic Macrophage in Ovarian Aging
## 从每个样本目录中读入表达矩阵文件并生成Seurt对象**********************************************
message(paste0(Sys.time(), "Single-Cell Atlas of Human Ovaries Reveals The Role Of The Pyroptotic Macrophage in Ovarian Aging"))
path <- "/root/wangje/Project/OveryArtical/Aritical2"
flist <- setNames(lapply(seq(421447, 421453), FUN = function(x) {
  message(paste0("HRS", x))
  RenameCells(CreateSeuratObject(Read10X(file.path(path, paste0("HRS", x), paste0("HRS", x), "outs", "filtered_feature_bc_matrix")),
    min.cell = 3, min.feature = 200
  ), add.cell.id = paste0("HRS", x))
}), paste0("HRS", seq(421447, 421453)))
## 合并Seurat对象列表生成一个Seurat对象*********************************************************
for (name in names(flist)) {
  flist[[name]][["sample"]] <- name
}
# scRNA.raw <- merge(flist[[1]],flist[-1])
# message(sprintf(c("基因为:%s","细胞为:%s"),dim(scRNA.raw)))

#**********************************************************************************************
# Spatiotemporal transcriptomic changes of human ovarian aging and the regulatory role of FOXP1
message(paste0(Sys.time(), "Spatiotemporal transcriptomic changes of human ovarian aging and the regulatory role of FOXP1"))
path <- "/root/wangje/Project/OveryArtical/GSE255690/scRNA/ExpressionData"
young <- setNames(lapply(file.path(path, "young", c("young1", "young3", "young4")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("young1", "young3", "young4"))

middle <- setNames(lapply(file.path(path, "middle", c("middle2", "middle3", "middle4")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("middle2", "middle3", "middle4"))

old <- setNames(lapply(file.path(path, "old", c("old1", "old2", "old3")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("old1", "old2", "old3"))
# 合并列表
flist1 <- c(young, middle, old)
# 添加sample列
for (name in names(flist1)) {
  flist1[[name]][["sample"]] <- name
  flist1[[name]] <- RenameCells(flist1[[name]], add.cell.id = name)
}

#**********************************************************************************************
# 自测数据
message(paste0(Sys.time(), "自测数据"))
inputDir <- "/root/wangje/Project/Data/rawData"
setwd(inputDir)
dirs <- list.files("./")
flist2 <- list()
for (i in dirs) {
  for (j in list.files(i)) {
    tryCatch(
      {
        message(paste0(Sys.time(), " ", i, " ", j))
        print(file.path(inputDir, i, j, "filtered_feature_bc_matrix/"))
        tmp <- Seurat::Read10X(file.path(inputDir, i, j, "filtered_feature_bc_matrix/"))
        tmp <- Seurat::CreateSeuratObject(tmp, min.cell = 3, min.feature = 200)
        tmp <- Seurat::RenameCells(tmp, add.cell.id = j)
        flist2[[paste0(i, "|", j)]] <- tmp
      },
      error = function(e) {
        message(e$message)
        return(NULL)
      }
    )
  }
}
for (i in names(flist2)) {
  flist2[[i]]$sample <- stringr::str_split_fixed(names(flist2), "\\|", n = 2)[, 2]
}
#***********************************************************************************************
# 合并数据
AllList <- c(flist, flist1, flist2)
scRNA <- merge(AllList[[1]], AllList[-1])
message(sprintf(c("基因为:%s", "细胞为:%s"), dim(scRNA)))
scRNA <- PercentageFeatureSet(scRNA, pattern = '^MT', col.name = 'percent.mt')
#***********************************************************************************************
setwd("/root/wangje/Project/OveryArtical")
qsave(scRNA, file = "./01_mergeRawSeurat.qs")
# 绘制小提琴图
scRNA$group <- ifelse(scRNA$sample %in% c("HRS421447", "HRS421448", "HRS421449", "HRS421450", "HRS421451", "HRS421452", "HRS421453"), "Artical1", ifelse(
  scRNA$sample %in% c("young1", "young3", "young4", "middle2", "middle3", "middle4", "old1", "old2", "old3"), "Artical2", "Self"
))
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"01_rawSeuratOnjectVlnplot.png"), height = 8, width = 16,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)

# 写出表格信息
message(paste0(Sys.time(), " 绘制统计结果... ..."))
df <- as.data.frame(scRNA@meta.data %>% dplyr::select(sample, group) %>% distinct())
rownames(df) <- 1:nrow(df)
colnames(df) <- c("Patient", "Group")
df$nCells <- sapply(unique(df$Patient), FUN = function(x) {
  ncol(scRNA[, scRNA$sample == x])
})
df$nGenes <- sapply(unique(df$Patient), FUN = function(x) {
  nrow(scRNA[, scRNA$sample == x])
})
df$`nCount_RNA(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(min(scRNA[, scRNA$sample == x]$nCount_RNA), "->", max(scRNA[, scRNA$sample == x]$nCount_RNA))
})
df$nCount_RNA_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$nCount_RNA)
})
df$`nFeature_RNA(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(min(scRNA[, scRNA$sample == x]$nFeature_RNA), "->", max(scRNA[, scRNA$sample == x]$nFeature_RNA))
})
df$nFeature_RNA_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$nFeature_RNA)
})
df$`percent.mt(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(round(min(scRNA[, scRNA$sample == x]$percent.mt), digits = 4), "->", round(max(scRNA[, scRNA$sample == x]$percent.mt), digits = 4))
})
df$percent.mt_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$percent.mt)
})

# 绘制表格
tbody.style <- tbody_style(
  color = "black",
  fill = c("#e8f3de", "#d3e8bb"), hjust = 1, x = 0.9
)
ggtexttable(df,
  rows = NULL,
  theme = ttheme(
    colnames.style = colnames_style(color = "white", fill = "#8cc257"),
    tbody.style = tbody.style
  )
) -> t_table
# 保存统计表格
ggsave(
  filename = file.path('./', "01_rawSeuratOnjectCount.png"), height = 8, width = 16,
  plot = t_table, limitsize = FALSE, bg = "white"
)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  数据粗过滤并使用DoubletFinder检测双细胞
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
filter_info <- scRNA$nCount_RNA <30000 & scRNA$nFeature_RNA > 200 & scRNA$nFeature_RNA < 7000 & scRNA$percent.mt < 20
scRNA.filter <- scRNA[,filter_info]
scRNA.filter$group <- ifelse(scRNA.filter$sample %in% c("HRS421447", "HRS421448", "HRS421449", "HRS421450", "HRS421451", "HRS421452", "HRS421453"), "Artical1", ifelse(
  scRNA$sample %in% c("young1", "young3", "young4", "middle2", "middle3", "middle4", "old1", "old2", "old3"), "Artical2", "Self"
))
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA.filter, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA.filter, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"02_rawFilterSeuratOnjectVlnplot.png"), height = 8, width = 16,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)
# scRNA.filter 33047 features across 227793 samples within 1 assay
qsave(scRNA.filter, file = "./02_rawFilterSeurat.qs")
#进行双细胞检测*************************************************************************************
rm(scRNA);gc()
Merge_obj <- NormalizeData(scRNA.filter, normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 3000)
Merge_obj <- ScaleData(Merge_obj, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
Merge_obj <- RunPCA(Merge_obj, npcs = 50, verbose = FALSE)

#去批次
library(harmony)
Merge_obj <- RunHarmony(Merge_obj,"sample")
Merge_obj <- RunUMAP(Merge_obj,  dims = 1:20, reduction = "harmony")
Merge_obj <- FindNeighbors(Merge_obj, reduction = "harmony",dims = 1:20) 
Merge_obj  <- FindClusters(object = Merge_obj , resolution = seq(from = 0.1, to = 1.0,  by = 0.1))

#单个分开，用来做DoubletFinder
sce_list <- SplitObject(Merge_obj, split.by = "sample")

detectDoublet <- function(
  obj, #seurat obj
  dims, #Pc数目
  estDubRate, #期望双细胞率，该值最好通过细胞在10X/Drop-Seq装置上的负载密度来估计
  ncores=1,#线程数
  SCTransform,#True or False
  Homotypic=F,#是否需要优化同源双细胞
  annotation){ #听说最好是celltype
  #use DoubletFinder packages
  require(DoubletFinder)#2.0.4
  
  #select pK
  sweep.res.list <- paramSweep(obj, PCs=dims, sct=SCTransform, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  
  
  #Doublet Proportion Estimate
  if(Homotypic==F){
    
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    
  }else{
    
    homotypic.prop <- modelHomotypic(obj@meta.data[,annotation])
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    nExp_poi  <- round(nExp_poi*(1-homotypic.prop))
  }
  
  # DoubletFinder:
  obj <- doubletFinder(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = SCTransform)
  
  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL
  
  return(obj)
}

# 判断预估的双细胞比例（10x的表格）
estDubRateFun <- function(srt){
    if(!inherits(srt,"Seurat")){stop("input file is not Seurat Object...")}
    rate <- tryCatch({
        if(dim(srt)[2] > 500 && dim(srt)[2] <=1000){
            tmp <- 0.004 + ((dim(srt)[2]-500) * round(0.004/500,digits = 7))
        }else if (dim(srt)[2] > 1000 && dim(srt)[2] <= 2000) {
            tmp <- 0.008 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 2000 && dim(srt)[2] <= 3000) {
            tmp <- 0.016 + ((dim(srt)[2]-1000) * round(0.007/1000,digits = 7))
        }else if (dim(srt)[2] > 3000 && dim(srt)[2] <= 4000) {
            tmp <- 0.023 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 4000 && dim(srt)[2] <= 5000) {
           tmp <- 0.031 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 5000 && dim(srt)[2] <= 6000) {
           tmp <- 0.039 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 6000 && dim(srt)[2] <= 7000) {
           tmp <- 0.046 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 7000 && dim(srt)[2] <= 8000) {
           tmp <- 0.054 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 8000 && dim(srt)[2] <= 9000) {
           tmp <- 0.061 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 9000 && dim(srt)[2] <= 10000) {
           tmp <- 0.069 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if(dim(srt)[2] >= 10000){
           tmp <- 0.076 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else {
           return(NULL)
        }
        tmp
    },error = function(e){
        message(e$error)
        return(NULL)
    })
    return(rate)
}

# 对每个样本进行分析
message(paste0(Sys.time(),' 开始进行双细胞检测...'))
sce_list <- lapply(sce_list, FUN = function(x){
    detectDoublet(x,dims = 1:30,estDubRate=estDubRateFun(x),
                                       ncores = 2, SCTransform=F, Homotypic=F, annotation="seurat_clusters")
})
message(paste0(Sys.time(),' 保存DoubletFinder结果文件...'))
saveRDS(sce_list, file = file.path(outputDir,'03_doubletFinderResult.rds'))

# 进行绘图
plist <- lapply(names(sce_list),FUN = function(x){
    Seurat::DimPlot(sce_list[[x]], group.by = "DF.classify") + labs(title = x)
})

png(file.path('./','03_doubletFinderResult.png'),res = 300, units = 'in', height = 4, width = 36)
print(patchwork::wrap_plots(plist, ncol = 9))
dev.off()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  使用R包SCP进行去批次和聚类
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 去除预测的双细胞
scRNA  <- merge(sce_list[[1]],sce_list[-1])
scRNA$new_group <- paste0(scRNA$sample,'_',scRNA$DF.classify)
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA.filter, group.by = "new_group", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA.filter, group.by = "new_group", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"03_DoubletFinderSeuratOnjectVlnplot.png"), height = 8, width = 26,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)
scRNA <- scRNA[,scRNA$DF.classify=="Singlet"]
library(magrittr)
scRNA %<>% NormalizeData()
scRNA %<>% CellCycleScoring(g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes,
                                  g1.features = cc.genes$g1.genes)
scRNA <- PercentageFeatureSet(scRNA, pattern = '^MT-', col.name = 'percent.mt') 
# harmony重新聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'Harmony',vars.to.regress=c('percent.mt','nCount_RNA','Phase'),cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'BBKNN',vars.to.regress=c('percent.mt','nCount_RNA','Phase'),cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'scVI',vars.to.regress=c('percent.mt','nCount_RNA','Phase'),cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  使用harmony同时去除两个批次(sample and group)
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(magrittr)
scRNA <- Seurat::CreateSeuratObject(counts = scRNA[["RNA"]]@counts, meta.data = scRNA@meta.data)
scRNA_harmony %<>% Seurat::NormalizeData()
scRNA_harmony %<>% CellCycleScoring(
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes,
  g1.features = cc.genes$g1.genes
) %<>% FindVariableFeatures() %<>% ScaleData(vars.to.regress = c("Phase", "nCount_RNA", "percent.mt")) %<>% RunPCA(verbose = FALSE)
system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = c("sample","group"))})
plot1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by = "sample", raster = TRUE)
plot2 <- ElbowPlot(scRNA_harmony, ndims = 50, reduction = "pca")
plotc <- plot1 + plot2
ggsave(filename = file.path('./', paste0("大群harmony", "ElbowPlot.png")), height = 4, width = 9, plot = plotc, bg = "white")
pc.num <- 1:40
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
#scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = pc.num)
scRNA_harmony <- FindClusters(scRNA_harmony, reduction = "harmony", resolution = seq(0.4,1,0.2))

# DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, raster = FALSE)
#DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
qsave(scRNA_harmony, file.path('./', paste0('001', "_大群Harmony.qs")))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  添加信息
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# singleR粗注释与
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
message(paste0('[',Sys.time(),']','进行singleR注释......'))
suppressMessages({
    library(Seurat)
    library(dplyr)
    library(patchwork)
    library(ggplot2)
    library(SingleR)
})

filePrefix <- "001AllCells"
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

scRNA <- Run_singleR(scRNA,cluster = scRNA$RNA_snn_res.1,resolution = "RNA_snn_res.1")

message(paste0('[',Sys.time(),']','开始进行绘图......'))
### DotPlot ###
marks <- function(){
    mark_list <- list(
    Immune = c("PTPRC"),
    NK = c("KLRC1", "FCGR3A", "NCAM1", "KLRD1", "GNLY"),
    "T cells" = c("CD3D", "CD3G", "CD2", "CD8A", "CD8B", "CD4"),
    Fibroblasts = c("COL1A1", "DCN", "LUM", "PDGFRA"),
    Myeloids = c("LYZ", "CD68", "TYROBP"),
    Epithelial = c("CD24", "KRT19", "EPCAM", "KRT18"),
    Bcells = c("CD79A", "CD19", "MS4A1"),
    Endothelial = c("CLDN5", "FLT1", "RAMP2", "CDH5"),
    Adipocytes = c('FABP4','PPARG','PLIN1'),
    Erythrocytes = c('HBB','HBA1','HBG1','HBD'),
    "Smooth muscle" = c("GJA4", "PLAC9", "CRYAB", "ACTA2", "PLN", "ADIRF", "MYH11", "CNN1", "MYL9"),
    Neutrophil = c("FCGR3B", "S100P", "CMTM2"),
    oocytes = c("TUBB8", "ZP3", "FIGLA"),
    Granulosal = c("GSTA1", "AMH", "HSD17B1", "DSP"),
    DC = c("LILRA4", "CXCR3", "IRF7"),
    Mast = c("CPA3", "TPSAB1", "TPSB2"))
    return(mark_list)
}
# 气泡图
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

raw_plot <- function(srt,resolution="RNA_snn_res.0.8",filePrefix){
    p1 <- DimPlot(srt, group.by = resolution, raster = TRUE,label = TRUE, repel = TRUE) + ggtitle(paste0('nCells: ',dim(srt)[2]))
    p2 <- DimPlot(srt, group.by = "singleR", raster = TRUE,label = TRUE, repel = TRUE)
    p3 <- DimPlot(srt, group.by = "sample", raster = FALSE,label = TRUE, repel = TRUE)
    p4 <- DimPlot(srt, group.by = "group", raster = FALSE,label = TRUE, repel = TRUE)
    ggsave(filename = paste0(filePrefix,"_singleR.png"),height = 4,width = 20,plot = p1+p2+p3+p4+patchwork::plot_layout(widths = c(1,1,1.5,1))+patchwork::plot_annotation(tag_levels = "A"))
    # 气泡图
    p5 <- plotBigDotPlot(srt, group.by = "singleR", marker = marks())
    p6 <- plotBigDotPlot(srt, group.by = resolution, marker = marks())
    ggsave(
    filename = paste0(filePrefix, "_singleR注释结果.png"), height = 14, width = 45, limitsize = F,
    plot = cowplot::plot_grid(plotlist = list(p5,p6), align = "hv", ncol = 2, nrow = 2, rel_widths = c(1, 1)), bg = "white")
}

raw_plot(srt = scRNA,filePrefix = filePrefix)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  细胞命名 进一步进行分群
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
scRNA@meta.data <- scRNA@meta.data %>%
  mutate(
    celltype = case_when(
      grepl('CD8\\+ T-cells|T cells|NK_cell|NK & T cells', singleR) ~ 'NK & T cells',
      grepl('Monocytes|Basophils|Macrophages|Neutrophils|Myeloid cells', singleR) ~ 'Myeloid cells',
      grepl('Mesangial cells|Epithelial_cells|Epithelial cells', singleR) ~ 'Epithelial cells',
      grepl('B-cells|B & Plasma cells', singleR) ~ 'B & Plasma cells',
      grepl('Erythroid cells|Erythrocytes', singleR) ~ 'Erythroid cells',
      grepl('Adipocytes|Endothelial_cells', singleR) ~ 'Endothelial cells',
      TRUE ~ singleR
    )
  )

scRNA$celltype <- ifelse(scRNA$RNA_snn_res.0.8 == 23,'oocytes',ifelse(scRNA$RNA_snn_res.0.8 == 26,'Granulosal',scRNA$celltype))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  每一个分群进行harmony去批次分群
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#### 对大群中注释的亚群数据进行进一步分析 ####
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(harmony)
    library(ggplot2)
    library(ggpubr)
    library(cowplot)
    library(SingleR)
    library(patchwork)
    library(qs)
})

# 读入数据
# Test_data <- "/root/wangje/Project/OveryArtical/001_大群Harmony.qs"
celltype_column <- "celltype"
message(paste0('[',Sys.time(),']',"读入数据....."))
# scRNA <- qs::qread(Test_data)
if(!(celltype_column %in% colnames(scRNA@meta.data))){
    print("细胞注释的列不在读入的Seurat对象中")
}
message(sprintf("存在细胞类型%s",paste0(unique(scRNA@meta.data[[celltype_column]]),collapse = "\t")))
celltypes <- unique(scRNA@meta.data[[celltype_column]])
scRNA$celltype <- scRNA@meta.data[[celltype_column]]

#### Cluster ####
analyze_subcluster <- function(srt, celltypes, vars.to.regress=c('Phase', 'nCount_RNA', 'percent.mt'),output_dir = ".") {
    if(length(celltypes) >1){
        stop('celltype的长度为1')
    }
    message(paste0('[', Sys.time(), ']', "分亚群进行分析....."))
    message(sprintf("分亚群进行分析.....%s", celltypes))
    
    scRNA_sub <- subset(srt, subset = celltype == celltypes)
    # message(sprintf("%s的数据维度为%s,%s", celltypes, paste0(c('基因', '细胞'),dim(scRNA_sub), collapse = ',')))
    print(dim(scRNA_sub))
    message(paste0('[', Sys.time(), ']', '进行harmomny数据分析'))
    
    scRNA_sub <- Seurat::CreateSeuratObject(counts = scRNA_sub[['RNA']]@counts, meta.data = scRNA_sub@meta.data)
    scRNA_harmony <- NormalizeData(scRNA_sub)
    scRNA_harmony <- FindVariableFeatures(scRNA_sub)
    scRNA_harmony <- ScaleData(scRNA_harmony, vars.to.regress = vars.to.regress)
    scRNA_harmony <- RunPCA(scRNA_harmony, verbose = FALSE)
    
    system.time({
        scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = c("sample","group"))
    })
    plot1 <- DimPlot(scRNA_harmony, reduction = "pca", group.by = "sample", raster = TRUE)
    plot2 <- ElbowPlot(scRNA_harmony, ndims = 50, reduction = "pca")
    plotc <- plot1 + plot2
    
    ggsave(filename = file.path(output_dir, paste0(celltypes, "ElbowPlot.png")), height = 4, width = 9, plot = plotc, bg = "white")
    
    pc.num <- 1:30
    scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = pc.num)
    #scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = pc.num)
    scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = pc.num)
    scRNA_harmony <- FindClusters(scRNA_harmony, reduction = "harmony", resolution = 0.8)
    
    DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, raster = FALSE)
    #DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
    qsave(scRNA_harmony, file.path(output_dir, paste0(celltypes, "_Harmony.qs")))
    return(scRNA_harmony)
}   

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

#### plot ####
marks <- function(){
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
    Mast = c("CPA3", "TPSAB1", "TPSB2"))
    return(mark_list)
}
# 气泡图
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

raw_plot <- function(srt,resolution="RNA_snn_res.0.8",filePrefix){
    p1 <- DimPlot(srt, group.by = resolution, raster = TRUE,label = TRUE, repel = TRUE) + ggtitle(paste0('nCells: ',dim(srt)[2]))
    p2 <- DimPlot(srt, group.by = "singleR", raster = TRUE,label = TRUE, repel = TRUE)
    p3 <- DimPlot(srt, group.by = "sample", raster = FALSE,label = TRUE, repel = TRUE)
    p4 <- DimPlot(srt, group.by = "group", raster = FALSE,label = TRUE, repel = TRUE)
    ggsave(filename = paste0(filePrefix,"_singleR.png"),height = 4,width = 20,plot = p1+p2+p3+p4+patchwork::plot_layout(widths = c(1,1,1.5,1))+patchwork::plot_annotation(tag_levels = "A"))
    # 气泡图
    p5 <- plotBigDotPlot(srt, group.by = "singleR", marker = marks())
    p6 <- plotBigDotPlot(srt, group.by = resolution, marker = marks())
    ggsave(
    filename = paste0(filePrefix, "_singleR注释结果.png"), height = 14, width = 30, limitsize = F,
    plot = cowplot::plot_grid(plotlist = list(p5,p6), align = "hv", ncol = 2, nrow = 2, rel_widths = c(1, 1)), bg = "white")
}


message(paste0("进行数据分群和SingleR粗注释......"))
result = NULL
output_dir = "."
for(cell in unique(celltypes)){
    print(paste0('[',Sys.time(),']',cell,' 使用harmony进行分群.....'))
    result <- analyze_subcluster(
        scRNA,celltypes = cell
    )
    # 分群结束，使用是singleR进行粗注释
    result = Run_singleR(data = result, cluster = result$RNA_snn_res.0.8, resolution = "RNA_snn_res.0.8")
    raw_plot(srt = result,filePrefix = cell)
    qsave(result, file.path(output_dir, paste0(cell, "_Harmony.qs")))
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  筛选出celltype中其他类型的细胞
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(ggplot2)
    library(Seurat)
    library(qs)
})
setwd('/root/wangje/Project/OveryArtical/New')
##### B cells ######
scRNA <- qread("B & Plasma cells_Harmony.qs") #33047 features across 1562 samples within 1 assay 
# B细胞中混有少量T细胞，DC细胞和成纤维细胞，其中成纤维细胞可能是双细胞，不建议后续使用
qsave(subset(scRNA,subset=RNA_snn_res.0.8==5),file = "B细胞中筛选DC细胞")
qsave(subset(scRNA,subset=RNA_snn_res.0.8==12),file = "B细胞中筛选T细胞")
qsave(subset(scRNA,subset=RNA_snn_res.0.8==7),file = "B细胞中筛选SMC细胞")
qsave(subset(scRNA,subset=RNA_snn_res.0.8==14),file = "B细胞中筛选Fibrobalsts细胞")
qsave(scRNA[,!scRNA$RNA_snn_res.0.8 %in% c(5,7,12,14)],file = "B细胞中筛选B细胞")

##### Myeloids ######
scRNA <- qread('Myeloid cells_Harmony.qs') # 33047 features across 13121 samples within 1 assay 
scRNA <- scRNA[,!scRNA$RNA_snn_res.0.8 %in% c(6,20,23,24)]
qsave(scRNA,file = "Myeloid cells_Harmony_sub.qs")

##### T cells ######
# 其中有部分SMC细胞可以去除
scRNA <- qread('NK & T cells_Harmony.qs') #33047 features across 20220 samples within 1 assay 
qsave(scRNA[,scRNA$RNA_snn_res.0.8==13],file = "NK & T cells_Harmony筛选SMC.qs") 
qsave(scRNA[,!scRNA$RNA_snn_res.0.8==13],file = "NK & T cells_Harmony筛选T细胞.qs") 

#### Fibroblasts ######
scRNA <- qread('/root/wangje/Project/OveryArtical/New/Fibroblasts_Harmony.qs') # 33047 features across 157935 samples within 1 assay
qsave(scRNA[,scRNA$RNA_snn_res.0.8 != 30],file = '/root/wangje/Project/OveryArtical/New/Fibroblasts_Harmony_sub.qs')

#### Endothelials #####
scRNA <- qread('/root/wangje/Project/OveryArtical/New/Endothelial cells_Harmony.qs') #33047 features across 25442 samples within 1 assay 
qsave(scRNA[,!scRNA$RNA_snn_res.0.8 %in% c(0,10,18)],file = '/root/wangje/Project/OveryArtical/New/Endothelial cells_Harmony_sub.qs')

#### Epithelials #####
scRNA <- qread('/root/wangje/Project/OveryArtical/New/Epithelial cells_Harmony.qs') # 33047 features across 9466samples within 1 assay
qsave(scRNA[,!scRNA$RNA_snn_res.0.8 %in% c(21,11,5,13)],file = '/root/wangje/Project/OveryArtical/New/Epithelial cells_Harmony_sub.qs')

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  筛选后的细胞再次进行聚类分群
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
library(Seurat)
library(dplyr)
library(qs)
library(SCP)


path <- '/root/wangje/Project/OveryArtical/New/'
getData <- function(filepath){
  tmp  <- qread(filepath)
  tmp <- Seurat::CreateSeuratObject(tmp@assays$RNA@counts,meta.data = tmp@meta.data)
  tmp
}

message("读入数据")
files <- list(
  # Myeloids = getData(paste0(path,'Myeloid cells_Harmony_sub.qs')),
  Tcells = getData(paste0(path,'NK & T cells_Harmony筛选T细胞.qs')),
  Endothelial = getData('/root/wangje/Project/OveryArtical/New/Endothelial cells_Harmony_sub.qs'),
  Epithelial = getData('/root/wangje/Project/OveryArtical/New/Epithelial cells_Harmony_sub.qs')
  # Fibroblasts = getData(paste0(path,'Fibroblasts_Harmony_sub.qs'))
)

files <- lapply(files,FUN = function(x){ 
  # 去除细胞数小于20的样本
  tmp  <- x
  tmp <- tmp[,tmp$sample %in% names(table(tmp$sample))[table(tmp$sample) >=20]]})

message("聚类分析") 
for(cell in names(files)){
  message(paste0(Sys.time(),'开始分析',cell)) 
  # scRNA <- SCP::Integration_SCP(files[[cell]], batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
  # qsave(scRNA, file = paste0(cell,"_scp.qs"))
  # BBKNN聚类
  scRNA <- SCP::Integration_SCP(files[[cell]], batch = "sample", integration_method = "BBKNN", cluster_resolution = seq(0.1, 1.5, 0.1))
  qsave(scRNA, file = paste0(cell,"_scp.qs"))
  # scVI聚类
  scRNA <- SCP::Integration_SCP(files[[cell]], batch = "sample", integration_method = "scVI", cluster_resolution = seq(0.1, 1.5, 0.1))
  qsave(scRNA, file = paste0(cell,"_scp.qs"))
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  进行绘图
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
