#install.packages("rio")

#加载包
library(tidyverse)
library(ggrepel)   # 用于标记
library(reshape2)
library(rio)

mutiVolcano = function(df,         
                       P = 0.05,   
                       FC = 1.414213562,
                       effect_size_threshold = 0.8,   # 新增：effect_size阈值
                       GroupName = c("Upregulated", "Downregulated", "Not Significant"),
                       pointColor = c("red", "blue", "gray"),
                       barFill = "#efefef",
                       tileLabel = "Label",
                       tileColor = RColorBrewer::brewer.pal(length(unique(df$Cluster)), "Set3")){
  # 修复FC参数与实际列名冲突问题
  FC_threshold = FC
  
  # 检查数据中是否有effect_size列
  if(!"effect_size" %in% colnames(df)) {
    stop("数据中缺少'effect_size'列。请确保输入数据包含effect_size值。")
  }
  
  # 数据分组（添加点尺寸和绘制顺序控制）
  dfSig = df %>%
    mutate(log2FC = log2(FC)) %>%
    filter(FC > FC_threshold | FC < 1/FC_threshold) %>%  # 修正筛选逻辑
    mutate(
      Group = case_when(
        log2FC > 0 & PValue < P & abs(effect_size) > effect_size_threshold ~ GroupName[[1]],
        log2FC < 0 & PValue < P & abs(effect_size) > effect_size_threshold ~ GroupName[[2]],
        TRUE ~ GroupName[[3]]
      ),
      PointShape = ifelse(VIP > 1, "VIP>1", "VIP<1"),
      PointSize = ifelse(Group == GroupName[[3]], 0.5, 1.5),  # 尺寸差异
      PlotOrder = ifelse(Group == GroupName[[3]], 1, 2)       # 绘制顺序控制
    ) %>%
    mutate(
      Group = factor(Group, levels = GroupName),
      Cluster = factor(Cluster, levels = unique(Cluster))
    ) %>%
    arrange(PlotOrder)  # 确保灰点先绘制
  
  # 背景柱形图数据
  global_min = -2
  global_max = 4.5
  dfBar = dfSig %>% group_by(Cluster) %>% summarise(min = global_min, max = global_max)
  
  # 散点图数据（添加抖动）
  dfJitter = dfSig %>% mutate(jitter = jitter(as.numeric(Cluster), factor = 2))
  
  # 绘图逻辑
  ggplot() +
    # 背景柱形
    geom_col(data = dfBar, aes(x = Cluster, y = max), fill = barFill) +
    geom_col(data = dfBar, aes(x = Cluster, y = min), fill = barFill) +
    
    # 点图层（关键修改）
    geom_point(
      data = dfJitter,
      aes(x = jitter, y = log2FC, color = Group, shape = PointShape, size = PointSize),
      stroke = 0.5
    ) +
    
    # 视觉映射系统
    scale_shape_manual(
      name = "VIP",
      values = c("VIP>1" = 17, "VIP=1" = 17, "VIP<1" = 16),
      labels = c("VIP>1" = "VIP ≥ 1", "VIP<1" = "VIP < 1")
    ) +
    scale_color_manual(values = setNames(pointColor, GroupName)) +
    scale_size_identity(guide = "none") +  # 尺寸不显示图例
    
    # 中间标签方块
    geom_tile(
      data = dfSig,
      aes(x = Cluster, y = 0, fill = Cluster), 
      color = "black",
      height = log2(FC_threshold) * 1.5,
      show.legend = FALSE
    ) +
    
    # 添加图例说明
    labs(
      x = "Clusters", 
      y = "log2FC",
      caption = paste0("Significance criteria: P < ", P, 
                       ", |log2FC| > ", round(log2(FC_threshold), 2),
                       ", |effect_size| > ", effect_size_threshold)
    ) +
    
    # 主题美化
    theme_classic() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      plot.caption = element_text(hjust = 0, size = 9, color = "gray50")
    ) +
    guides(
      color = guide_legend(order = 1, title = "Expression"),
      shape = guide_legend(order = 2, title = "VIP Status")
    )
}
#加载数据
df <- read.csv("3soil_diff_mutiVolcano.csv", header=T)
pdf('3soil_diff_mutiVolcano.pdf',width = 10,height = 6)
mutiVolcano(
  df = df,    # 绘图数据
  P = 0.05,   # P值卡值
  FC = 1.414213562,   # FC卡值
  effect_size_threshold = 0.8,  # effect_size阈值
  GroupName = c("Up-regulated", "Down-regulated", "Not Significant"),
  pointColor = c("red", "blue", "gray"),
  barFill = "#efefef",
  tileLabel = "Label",
  tileColor = RColorBrewer::brewer.pal(length(unique(df$Cluster)), "Set3")
)
dev.off()

