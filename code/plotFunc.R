
suppressMessages({
    library(Seurat)
    library(dplyr)
    library(patchwork)
    library(dplyr)
    library(cowplot)
    library(tidydr)
    library(RColorBrewer)
    library(scCustomize)
    library(ggpubr)
    library(ggplot2)
    library(grid)
    library(gridExtra)
    library(BiocParallel)
    library(ComplexHeatmap)
    library(pheatmap)
    library(stringr)
    library(SeuratObject)
})

#坐标轴绘制，修改SCP中的theme_blank函数
theme_blank1 <- function(
    add_coord = TRUE, xlen_npc = 0.25, ylen_npc = 0.25,
    xlab = "", ylab = "", lab_size = 12, ...) {
  if (isTRUE(add_coord) && isTRUE(xlab != "")) {
    x_space <- lab_size + 2
  } else {
    x_space <- 0
  }
  if (isTRUE(add_coord) && isTRUE(ylab != "")) {
    y_space <- lab_size + 2
  } else {
    y_space <- 0
  }
  args1 <- list(
    panel.border = element_blank(), panel.grid = element_blank(),
    axis.title = element_blank(), axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text = element_blank(),
    legend.background = element_blank(), legend.box.margin = margin(
      0,
      0, 0, 0
    ), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(0,
      0, x_space, y_space,
      unit = "points"
    ), complete = FALSE
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(what = theme, args = args)
  if (isTRUE(add_coord)) {
    g <- grobTree(gList(linesGrob(x = unit(
      c(0, xlen_npc),
      "npc"
    ), y = unit(c(0, 0), "npc"), arrow = arrow(angle = 25, length = unit(
      0.045,
      "npc"
    ), type = "open"), gp = gpar(lwd = 3,fill = "black")), textGrob(
      label = xlab,
      x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3,
      hjust = 0, gp = gpar(fontsize = lab_size)
    ), linesGrob(x = unit(c(
      0,
      0
    ), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(angle = 25, length = unit(
      0.045,
      "npc"
    ), type = "open"), gp = gpar(lwd = 3,col="black")), textGrob(
      label = ylab,
      x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3,
      hjust = 0, rot = 90, gp = gpar(fontsize = lab_size)
    )))
    return(list(list(annotation_custom(g)), list(theme_scp() +
      out), list(coord_cartesian(clip = "off"))))
  } else {
    return(list(list(theme_scp() + out)))
  }
}
## 1 绘制DimPlot
#修改tidydr包中的theme_dr函数
theme_dr <- function (xlength = 0.3, ylength = 0.3, arrow = grid::arrow(length = unit(0.15,
    "inches"), type = "closed"))
{
    theme_minimal() %+replace% theme_noaxis(axis.line.x.bottom = element_line2(id = 1,
        xlength = xlength, arrow = arrow), axis.line.y.left = element_line2(id = 2,
        ylength = ylength, arrow = arrow), axis.title = element_text(hjust = 0.0))
}

theme_set <- theme_dr(arrow = grid::arrow(angle = 20, length = unit(0.15, "inches"), type = "closed"),xlength = 0.24, ylength = 0.24) +
    theme(plot.title = element_blank(), 
    panel.grid = element_blank(),
    legend.text = element_text(size = 12,color="black"),
    axis.title = element_text(size = 14,colour = "black")) 

## 2 Add information
library(magrittr)
# scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data),'_[A|T|G|C].*',n=2)[,1]
scRNA@meta.data %<>% dplyr::mutate(
    group1 = dplyr::case_when(
        sample %in% c('young1','young3','young4',"HRS421451","HRS421452",'na1111','HRS421447') ~'18-29',
        sample %in% c('na1122','HRS421453','na-426','na-1010') ~ '32-35',
        sample %in% c('middle2','middle3','middle4','HRS421448') ~ '37-39',
        sample %in% c('na1128','HRS421449','HRS421450','na1117') ~ '39-44',
        sample %in% c('old1','old2','old3') ~ '47-49',
        sample %in% c('na1129','NA_728','na-412') ~ '55-60'
        ),
    group1 = factor(group1,levels = c('18-29','32-35','37-39','39-44','47-49','55-60'))
)


