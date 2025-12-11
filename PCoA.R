# =============================
# Beta多样性分析 - 完整统计分析
# =============================

# 数据载入
rm(list=ls())

# 加载必要的包
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggthemes", quietly = TRUE)) install.packages("ggthemes")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")# 从GitHub安装ggalt
  devtools::install_github("hrbrmstr/ggalt")

library(vegan)
library(ggplot2)
library(ggthemes)
library(ggalt)
library(ggpubr)
library(patchwork)
library(pairwiseAdonis)

# =============================
# 数据准备
# =============================
# 加载物种丰度数据
spe <- read.csv("map00910_percentages.txt", header = TRUE, row.names = 1, sep = "\t")
totu <- t(spe)
otu <- as.data.frame(totu)

# 加载分组信息
group <- read.csv('group.csv', row.names = 1)
group_vector <- group$group  # 假设分组信息在group列

# =============================
# 1. 计算Bray-Curtis距离
# =============================
distance <- vegdist(otu, method = 'bray')

# =============================
# 2. PERMANOVA分析
# =============================
# 整体PERMANOVA
set.seed(123)  # 设置随机种子确保结果可重复
permanova_result <- adonis2(distance ~ group_vector, 
                            permutations = 999, 
                            method = "bray")

# 输出整体PERMANOVA结果
cat("\n========== 整体PERMANOVA分析结果 ==========\n")
print(permanova_result)

# 两两比较PERMANOVA
pairwise_permanova <- pairwise.adonis(
  x = otu, 
  factors = group_vector,
  sim.function = "vegdist",
  sim.method = "bray",
  p.adjust.m = "BH",  # 使用Benjamini-Hochberg FDR校正
  reduce = NULL,
  perm = 999
)

# 输出两两比较结果
cat("\n========== 两两比较PERMANOVA分析结果 ==========\n")
print(pairwise_permanova)

# =============================
# 3. Beta-dispersion分析（同质性检验）
# =============================
# 计算各组的多变量离散度
beta_disp <- betadisper(distance, group_vector)

# 置换检验离散度的同质性
disp_permutest <- permutest(beta_disp, 
                            pairwise = TRUE, 
                            permutations = 999)

# 输出离散度检验结果
cat("\n========== Beta-dispersion分析结果 ==========\n")
cat("整体离散度检验:\n")
print(disp_permutest$tab)

cat("\n两两离散度比较:\n")
print(disp_permutest$pairwise)

# 计算各组离散度统计
disp_stats <- data.frame(
  Group = levels(factor(group_vector)),
  Average_Distance_to_Centroid = beta_disp$distances,
  Group_Size = table(group_vector)
)

# =============================
# 4. 保存统计分析结果
# =============================
# 创建结果汇总列表
stat_summary <- list(
  overall_permanova = permanova_result,
  pairwise_permanova = pairwise_permanova,
  beta_dispersion = disp_permutest,
  dispersion_statistics = disp_stats
)

# 保存为R数据文件
saveRDS(stat_summary, file = "beta_diversity_statistics.rds")

# 保存为文本文件
sink("beta_diversity_statistics.txt")
cat("BETA DIVERSITY STATISTICAL ANALYSIS\n")
cat("====================================\n\n")

cat("1. OVERALL PERMANOVA\n")
cat("------------------------------------\n")
print(permanova_result)
cat("\n\n")

cat("2. PAIRWISE PERMANOVA (FDR-adjusted)\n")
cat("------------------------------------\n")
print(pairwise_permanova)
cat("\n\n")

cat("3. BETA-DISPERSION ANALYSIS\n")
cat("------------------------------------\n")
cat("Overall homogeneity test:\n")
print(disp_permutest$tab)
cat("\nPairwise dispersion comparisons:\n")
print(disp_permutest$pairwise)
sink()

# =============================
# 5. PCoA可视化
# =============================
# 执行PCoA
pcoa <- cmdscale(distance, k = (nrow(otu) - 1), eig = TRUE)

# 提取前两个主坐标
plot_data <- data.frame(PCoA1 = pcoa$point[,1], 
                        PCoA2 = pcoa$point[,2],
                        Group = group_vector,
                        ID = rownames(otu))

# 计算解释率
eig <- pcoa$eig
sum_eig <- sum(eig[eig > 0])  # 只考虑正特征值
eig_percent <- round(eig[1:2]/sum_eig * 100, 1)

# 创建统计信息字符串
stat_info <- paste0(
  "Overall PERMANOVA: F = ", round(permanova_result$F[1], 2), 
  ", R² = ", round(permanova_result$R2[1], 3), 
  ", p = ", round(permanova_result$`Pr(>F)`[1], 4),
  "\nBeta-dispersion: F = ", round(disp_permutest$tab$F[1], 2),
  ", p = ", round(disp_permutest$tab$`Pr(>F)`[1], 4)
)

# 创建主图
p <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2, color = Group, group = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2, size = 0.8, alpha = 0.5) +
  scale_fill_manual(values = c("#A5C496","#C7988C","#8891DB")) +  # 填充颜色
  scale_color_manual(values = c("#A5C496","#C7988C","#8891DB")) +  # 点颜色
  labs(
    x = paste0("PCoA 1 (", eig_percent[1], "%)"),
    y = paste0("PCoA 2 (", eig_percent[2], "%)"),
    title = "Principal Coordinates Analysis (PCoA)",
    subtitle = stat_info
  ) +
  geom_encircle(aes(fill=Group), alpha = 0.1, show.legend = F) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold")
  ) +
  coord_fixed(ratio = 1)

# =============================
# 6. 创建结果表格
# =============================
# 准备PERMANOVA两两比较表格
pairwise_table <- pairwise_permanova[, c("pairs", "F.Model", "R2", "p.value", "p.adjusted")]
names(pairwise_table) <- c("Comparison", "F", "R²", "p-value", "FDR-adjusted p")

# 创建表格图
tab <- ggtexttable(pairwise_table, 
                   rows = NULL,
                   theme = ttheme(
                     colnames.style = colnames_style(color = "white", fill = "#2E86C1"),
                     tbody.style = tbody_style(color = "black", fill = c("#F2F3F4", "#FFFFFF"))
                   )) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2) %>%
  tab_add_hline(at.row = nrow(pairwise_table) + 1, row.side = "bottom", linewidth = 2) %>%
  tab_add_title(text = "Pairwise PERMANOVA Comparisons (BH FDR-adjusted)", 
                face = "bold", size = 12)

# =============================
# 7. 组合图形
# =============================
# 使用patchwork组合图形
combined_plot <- p / tab + 
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = 'Beta Diversity Analysis: Microbial Community Structure',
    subtitle = paste('Based on Bray-Curtis dissimilarity |', 
                     'Total samples:', nrow(otu), '| Groups:', length(unique(group_vector))),
    caption = 'Statistical tests: PERMANOVA (999 permutations) with beta-dispersion testing',
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.caption = element_text(size = 10, hjust = 1, face = "italic")
    )
  )

ggsave("PCoA_complete_analysis.pdf", 
       plot = combined_plot,
       width = 12, 
       height = 14, 
       dpi = 300)

ggsave("PCoA_complete_analysis.png", 
       plot = combined_plot,
       width = 12, 
       height = 14, 
       dpi = 300, 
       bg = "white")

# =============================
# 8. 生成摘要报告
# =============================
cat("\n========== 分析完成 ==========\n")
cat("生成的文件:\n")
cat("1. beta_diversity_statistics.rds - R数据文件（完整结果）\n")
cat("2. beta_diversity_statistics.txt - 文本摘要文件\n")
cat("3. PCoA_complete_analysis.pdf - 完整分析图形（PDF）\n")
cat("4. PCoA_complete_analysis.png - 完整分析图形（PNG）\n")

cat("\n统计摘要:\n")
cat(sprintf("• 样本数量: %d\n", nrow(otu)))
cat(sprintf("• 组别数量: %d\n", length(unique(group_vector))))
cat(sprintf("• 整体PERMANOVA R²: %.3f (p = %.4f)\n", 
            permanova_result$R2[1], permanova_result$`Pr(>F)`[1]))
cat(sprintf("• Beta-dispersion检验p值: %.4f\n", disp_permutest$tab$`Pr(>F)`[1]))

# 显示图形
print(combined_plot)

