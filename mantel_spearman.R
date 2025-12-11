# =============================
# 加载必要包
# =============================
# 安装所需包（如果尚未安装）
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("linkET", quietly = TRUE))
  devtools::install_github("Hy4m/linkET", force = TRUE)

library(linkET)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 设置随机种子以确保可重复性
set.seed(123)

# =============================
# 数据加载与预处理
# =============================
# 读取数据
varechem <- read.csv("Properties.txt", header = TRUE, row.names = 1, sep = "\t")
varespec <- read.csv("k_Genus.txt", header = TRUE, row.names = 1, sep = "\t")

# 检查数据维度
cat("土壤理化性质数据维度：", dim(varechem), "\n")
cat("生物数据维度：", dim(varespec), "\n")

# 标准化土壤理化数据（Z-score标准化）
varechem_scaled <- scale(varechem)

# =============================
# Mantel检验分析
# =============================
# 执行Mantel检验 - 使用正确的参数格式
mantel_results <- mantel_test(
  varespec, 
  varechem_scaled,
  spec_select = list(
    "Nitrogen Metabolism Genes" = 1:53,
    "higher relative abundance Genera" = 54:73
  )
) 

# 添加统计信息列
mantel_processed <- mantel_results %>%
  mutate(
    # 计算95%置信区间
    CI_lower = r - 1.96 * sqrt((1 - r^2) / (nrow(varechem) - 2)),
    CI_upper = r + 1.96 * sqrt((1 - r^2) / (nrow(varechem) - 2)),
    # 对p值进行FDR校正
    q_value = p.adjust(p, method = "BH"),
    # 为可视化创建分类变量
    r_category = cut(r, 
                     breaks = c(-Inf, 0.3, 0.5, 0.7, Inf),
                     labels = c("Weak (<0.3)", "Moderate (0.3-0.5)", 
                                "Strong (0.5-0.7)", "Very Strong (>0.7)")),
    p_category = cut(q_value,  # 使用FDR校正后的q值
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "ns"))
  )

# 输出Mantel检验结果
cat("\n========== Mantel检验结果 ==========\n")
print(mantel_processed)

# 保存详细结果
write.csv(mantel_processed, "mantel_test_results.csv", row.names = FALSE)

# =============================
# 土壤理化性质间的相关性分析
# =============================
# 计算Spearman相关系数
# 使用 correlate() 函数并指定 is_corr = TRUE
soil_cor <- correlate(varechem_scaled, method = "spearman", is_corr = TRUE)

# 检查计算是否成功
if (!is.null(soil_cor)) {
  cat("相关性计算成功\n")
} else {
  cat("相关性计算失败，尝试替代方法...\n")
  # 如果 correlate() 失败，尝试其他方法
  soil_cor_matrix <- cor(varechem_scaled, method = "spearman")
  # 创建一个简单的相关性对象
  soil_cor <- list(r = soil_cor_matrix)
}

# =============================
# 整合可视化 - 修复版本
# =============================
# 首先创建一个相关性热图，设置 is_corr = TRUE
p1 <- qcorrplot(soil_cor$r, 
                type = "lower", 
                diag = FALSE,
                is_corr = TRUE) +  # 明确指定这是相关性矩阵
  # 添加相关系数方块
  geom_square() +
  # 添加显著性标记 - 简化版本
  geom_mark(sig_level = c(0.001, 0.01, 0.05),
            sig_thres = 0.05,
            size = 2.5,
            colour = "white") +
  # 添加Mantel检验连接线 - 修复数据引用
  geom_couple(data = mantel_processed,
              aes(colour = p_category, 
                  size = r_category),
              curvature = nice_curvature()) +
  # 颜色和大小标度
  scale_fill_gradientn(
    name = "Spearman's ρ",
    colours = rev(brewer.pal(11, "RdBu")),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  scale_colour_manual(
    name = "FDR-adjusted q-value",
    values = c("***" = "#2E8B57", "**" = "#7CCD7C", "*" = "#9BCD9B", "ns" = "gray70"),
    labels = c("***" = "< 0.001", "**" = "0.001-0.01", "*" = "0.01-0.05", "ns" = "≥ 0.05")
  ) +
  scale_size_manual(
    name = "Mantel's r",
    values = c("Weak (<0.3)" = 0.5, "Moderate (0.3-0.5)" = 1, 
               "Strong (0.5-0.7)" = 1.5, "Very Strong (>0.7)" = 2)
  ) +
  # 图例调整
  guides(
    fill = guide_colorbar(order = 1),
    colour = guide_legend(order = 2, 
                          override.aes = list(size = 2)),
    size = guide_legend(order = 3,
                        override.aes = list(colour = "black"))
  ) +
  # 主题设置
  theme_minimal() +
  theme(
    plot.margin = margin(20, 40, 20, 20),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 保存图形
ggsave(
  plot = p1,
  filename = 'mantel_correlation_plot.pdf',
  width = 14,
  height = 10,
  device = cairo_pdf
)

cat("\n========== 分析完成 ==========\n")
cat("1. Mantel检验结果已保存至: mantel_test_results.csv\n")
cat("2. 整合可视化图已保存至: mantel_correlation_plot.pdf\n")

# =============================
# 生成统计摘要
# =============================
cat("\n========== 统计摘要 ==========\n")

# Mantel检验统计摘要
cat("\nMantel检验统计摘要:\n")
summary_stats <- mantel_processed %>%
  group_by(spec) %>%  # 注意：列名可能是 "spec" 而不是 "spec_group"
  summarise(
    n_tests = n(),
    n_significant = sum(q_value < 0.05, na.rm = TRUE),
    mean_r = mean(r, na.rm = TRUE),
    sd_r = sd(r, na.rm = TRUE),
    min_r = min(r, na.rm = TRUE),
    max_r = max(r, na.rm = TRUE)
  )
print(summary_stats)

# 土壤相关性统计摘要 - 简化版本
if (exists("soil_cor_matrix")) {
  cat("\n土壤性质相关性统计摘要:\n")
  # 将相关性矩阵转换为长格式
  soil_cor_long <- as.data.frame(as.table(soil_cor_matrix))
  names(soil_cor_long) <- c("Var1", "Var2", "rho")
  
  # 移除对角线和对角线上方的重复项
  soil_cor_long <- soil_cor_long[as.numeric(soil_cor_long$Var1) < as.numeric(soil_cor_long$Var2), ]
  
  soil_summary <- soil_cor_long %>%
    summarise(
      n_correlations = n(),
      mean_rho = mean(rho, na.rm = TRUE),
      sd_rho = sd(rho, na.rm = TRUE)
    )
  print(soil_summary)
}

