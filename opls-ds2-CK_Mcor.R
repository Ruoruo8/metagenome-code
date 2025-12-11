# =============================
# 数据载入与包加载
# =============================
rm(list=ls())

library(ropls)
library(ggplot2)
library(ggsci)
library(Cairo)
library(tidyverse)
library(extrafont)
loadfonts()
library(mixOmics)
library(pROC)
library(caret)
library(readr)
library(pheatmap)
library(effectsize)

set.seed(123)  # 全局随机种子

# ============================================
# 代谢组数据预处理脚本
# 适用于：PCA_KB_JZL.csv
# 第一列为代谢物编号，第一行为样本名称
# ============================================
# 3. 读取数据
cat("正在读取数据文件...\n")
raw_data <- read.csv("PCA_JZL_KB.csv", row.names = 1, check.names = FALSE)
cat("数据维度：", nrow(raw_data), "个代谢物 ×", ncol(raw_data), "个样本\n")

# 4. 数据质量检查
cat("\n=== 数据质量检查 ===\n")
# 4.1 检查缺失值
missing_percent <- sum(is.na(raw_data)) / (nrow(raw_data) * ncol(raw_data)) * 100
cat("1. 缺失值比例：", round(missing_percent, 2), "%\n")

# 4.2 检查零值比例
zero_percent <- sum(raw_data == 0, na.rm = TRUE) / (nrow(raw_data) * ncol(raw_data)) * 100
cat("2. 零值比例：", round(zero_percent, 2), "%\n")

# 4.3 查看数据范围
cat("3. 原始数据范围：", min(raw_data, na.rm = TRUE), "到", max(raw_data, na.rm = TRUE), "\n")

# 6. 数据转换（对数转换）
cat("\n=== 数据转换 ===\n")
# 检查是否有零或负值，如果有需要加一个小的常数
if(any(raw_data <= 0)) {
  constant <- abs(min(raw_data[raw_data <= 0])) + 1
  cat("数据包含非正值，添加常数：", constant, "\n")
  raw_data <- raw_data + constant
}

# 对数转换（以2为底）
data_log <- log2(raw_data + 1)  # +1 防止对0取对数
cat("完成对数转换 (log2(x+1))\n")

# 7. 样本间标准化（中位数标准化）
cat("\n=== 样本间标准化 ===\n")
# 计算每个样本的中位数
sample_medians <- apply(data_log, 2, median)
# 计算全局中位数
global_median <- median(sample_medians)
cat("样本中位数范围：", round(min(sample_medians), 2), "到", round(max(sample_medians), 2), "\n")
cat("全局中位数：", round(global_median, 2), "\n")

# 进行中位数标准化
data_normalized <- data_log
for(sample in colnames(data_log)) {
  scaling_factor <- global_median / sample_medians[sample]
  data_normalized[, sample] <- data_log[, sample] * scaling_factor
}
cat("完成中位数标准化\n")

# 8. 代谢物缩放（Pareto Scaling - 代谢组学最常用）
cat("\n=== Pareto缩放 ===\n")
# 计算每个代谢物的标准差
feature_sd <- apply(data_normalized, 1, sd)
# 进行Pareto缩放：除以标准差的平方根
data_scaled <- data_normalized
for(i in 1:nrow(data_normalized)) {
  data_scaled[i, ] <- data_normalized[i, ] / sqrt(feature_sd[i])
}
cat("完成Pareto缩放（除以各自标准差的平方根）\n")

# 9. 均值中心化（通常作为最后一步）
cat("\n=== 均值中心化 ===\n")
# 对每个代谢物进行均值中心化（按行）
row_means <- apply(data_scaled, 1, mean)
data_centered <- data_scaled - row_means
cat("完成均值中心化\n")

# =============================
# 1. 数据载入
# =============================
otu <- t(data_centered)  # 转置成样本行，变量列
group <- read.csv("group_JZL_KB.csv", row.names = 1)


# =============================
# 2. 无监督 PCA 分析
# =============================
pca_result <- prcomp(otu, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x[,1:2])
pca_scores$group <- group$group
pca_scores$samples <- rownames(pca_scores)
pca_var <- summary(pca_result)$importance[2,1:2] * 100

p_pca <- ggplot(pca_scores, aes(x=PC1, y=PC2, color=group)) +
  geom_point(size=3) +
  stat_ellipse(level=0.95, linetype=2, size=0.5) +
  labs(x=paste0("PC1 (", round(pca_var[1],2), "%)"),
       y=paste0("PC2 (", round(pca_var[2],2), "%)"),
       title="Unsupervised PCA") +
  theme_bw() +
  theme(panel.grid=element_blank(), legend.position="bottom")
ggsave(filename = "PCA_JZL_CK.pdf", plot = p_pca, width = 6, height = 6, device = cairo_pdf)

# =============================
# 4. OPLS-DA 模型
# =============================
df1_oplsda <- opls(
  x = otu, 
  y = group$group, 
  predI = 1, 
  orthoI = NA, 
)

# 模型统计
model_stats <- data.frame(
  Component = "p1",
  R2X = df1_oplsda@modelDF["p1", "R2X"],
  R2X_cum = df1_oplsda@modelDF["p1", "R2X(cum)"],
  R2Y = df1_oplsda@modelDF["p1", "R2Y"],
  R2Y_cum = df1_oplsda@modelDF["p1", "R2Y(cum)"],
  Q2 = df1_oplsda@modelDF["p1", "Q2"],
  Q2_cum = df1_oplsda@modelDF["p1", "Q2(cum)"]
)
print(model_stats)

# =============================
# 5. 置换检验 (Permutation Test)
# =============================
perm_result <- opls(
  x = otu,
  y = group$group,
  predI = 1,
  orthoI = NA,
  permI = 1000,
)

perm_result

# ropls 自带置换图
plot(perm_result, typeVc = "permutation")
#OPLS-DA_JZL_CK_permutation.pdf
plot(perm_result, typeVc = "summary")
#OPLS-DA_JZL_CK_summary.pdf

# 保存置换结果
perm_stats <- data.frame(
  perm_R2Y = perm_result@summaryDF$R2Y,
  perm_Q2 = perm_result@summaryDF$Q2,
  p_R2Y = perm_result@summaryDF$pR2Y,
  p_Q2 = perm_result@summaryDF$pQ2
)
print(perm_stats)

# =============================
# 6. 重复交叉验证 (10x7-fold)
# =============================
repeated_cv_results <- data.frame(
  R2Y = numeric(10),
  Q2 = numeric(10),
  classification_error = numeric(10)
)

for (rep in 1:10) {
  train_idx <- createDataPartition(group$group, p=0.7, list=FALSE)
  train_data <- otu[train_idx,]
  train_group <- group$group[train_idx]
  test_data <- otu[-train_idx,]
  test_group <- group$group[-train_idx]
  
  cv_model <- opls(train_data, train_group, predI=1, orthoI=NA, crossvalI=7)
  
  # 预测
  pred_values <- predict(cv_model, newdata=test_data)
  
  # 分类
  pred_classes <- ifelse(pred_values > 0.5, levels(factor(group$group))[2], levels(factor(group$group))[1])
  
  confusion_mat <- confusionMatrix(
    factor(pred_classes, levels=levels(factor(group$group))),
    factor(test_group, levels=levels(factor(group$group)))
  )
  
  repeated_cv_results$R2Y[rep] <- cv_model@summaryDF$R2Y
  repeated_cv_results$Q2[rep] <- cv_model@summaryDF$Q2
  repeated_cv_results$p_R2Y[rep] <- cv_model@summaryDF$pR2Y
  repeated_cv_results$p_Q2[rep] <- cv_model@summaryDF$pQ2
  repeated_cv_results$classification_error[rep] <- 1 - confusion_mat$overall["Accuracy"]
}

cv_summary <- data.frame(
  metric = c("R2Y", "Q2", "p_R2Y", "p_Q2", "Classification_Error"),
  mean = c(mean(repeated_cv_results$R2Y, na.rm=TRUE),
           mean(repeated_cv_results$Q2, na.rm=TRUE),
           mean(repeated_cv_results$p_R2Y, na.rm=TRUE),
           mean(repeated_cv_results$p_Q2, na.rm=TRUE),
           mean(repeated_cv_results$classification_error, na.rm=TRUE)),
  sd = c(sd(repeated_cv_results$R2Y, na.rm=TRUE),
         sd(repeated_cv_results$Q2, na.rm=TRUE),
         sd(repeated_cv_results$p_R2Y, na.rm=TRUE),
         sd(repeated_cv_results$p_Q2, na.rm=TRUE),
         sd(repeated_cv_results$classification_error, na.rm=TRUE))
)
print(cv_summary)

# =============================
# 7. 绘制 OPLS-DA Score Plot + 置信椭圆 + 重复CV统计
# =============================
data <- as.data.frame(df1_oplsda@scoreMN)
if (!is.null(df1_oplsda@orthoScoreMN)) {
  o1 <- df1_oplsda@orthoScoreMN[,1]
  data$o1 <- o1
} else {
  data$o1 <- 0
}
data$group <- group$group
data$samples <- rownames(data)
x_lab <- df1_oplsda@modelDF[1, "R2X"] * 100
col <- c("#C7988C", "#8891DB")

p_oplsda <- ggplot(data, aes(x=p1, y=o1, color=group)) +
  geom_point(size=3) +
  stat_ellipse(level=0.95, linetype=2, size=0.5, aes(fill=group), alpha=0.2) +
  geom_vline(xintercept=0, lty="dashed", color="red") +
  geom_hline(yintercept=0, lty="dashed", color="red") +
  labs(
    x = paste0("P1 (", round(x_lab,2), "%)"),
    y = "to1",
    title = paste0(
      "OPLS-DA Score Plot\n",
      "Repeated CV: R²Y = ", round(cv_summary$mean[1],3), " ± ", round(cv_summary$sd[1],3),
      ", Q² = ", round(cv_summary$mean[2],3), " ± ", round(cv_summary$sd[2],3),
      ", Classification Error = ", round(cv_summary$mean[3],3), " ± ", round(cv_summary$sd[3],3)
    )
  ) +
  scale_color_manual(values=col) +
  scale_fill_manual(values=col) +
  theme_bw() +
  theme(panel.grid=element_blank(), legend.position="bottom")

ggsave("OPLS-DA_JZL_CK_CV.pdf", p_oplsda, width=6, height=6, device=cairo_pdf)

# =============================
# 8. 绘制重复CV指标柱状图
# =============================
cv_long <- repeated_cv_results %>%
  mutate(rep=1:n()) %>%
  pivot_longer(cols=c("R2Y", "Q2", "classification_error"), 
               names_to="metric", values_to="value")

p_cv <- ggplot(cv_long, aes(x=factor(rep), y=value, fill=metric)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~metric, scales="free_y") +
  labs(x="Repeat", y="Value", title="Repeated 7-fold CV Results") +
  scale_fill_manual(values=c("R2Y"="#56B4E9", "Q2"="#8891DB", "classification_error"="#C7988C")) +
  theme_bw() +
  theme(panel.grid=element_blank(), legend.position="none")

ggsave("Repeated_CV_Q2_Error.pdf", p_cv, width=6, height=4, device=cairo_pdf)

# =============================
# 9. 提取VIP + FDR整合
# =============================
# 提取VIP值
data_VIP <- df1_oplsda@vipVn


# 注意：OPLS-DA本身不提供p值，需通过统计检验计算（如t检验或Wilcoxon）
# 假设数据已转置为otu矩阵（行=样本，列=变量）
p_values <- apply(t(raw_data), 2, function(x) {
  t.test(x[1:5],x[6:10], paired = TRUE)$p.value  # 关键修改点：调用t.test并添加paired参数
})

fdr <- p.adjust(p_values, method = "fdr")

effect_sizes <- apply(t(raw_data), 2, function(x) {
  cohens_d(x ~ group$group)$Cohens_d
})

# 计算log2 Fold Change
calculate_log2fc <- function(data_matrix, group_vector, pseudo_count = 1e-10) {
  # 添加伪计数避免log(0)
  data_adj <- data_matrix + pseudo_count
  
  # 计算组均值
  group_means <- aggregate(t(data_adj), 
                           list(group = group_vector), 
                           mean)
  
  # 转置回来
  group_means_t <- as.data.frame(t(group_means[, -1]))
  colnames(group_means_t) <- group_means$group
  
  # 计算log2FC（第二组/第一组）
  groups <- colnames(group_means_t)
  log2fc <- log2(group_means_t[, groups[2]] / group_means_t[, groups[1]])
  
  return(log2fc)
}

# 应用到您的数据
log2fc_values <- calculate_log2fc(raw_data, group$group)

# 合并VIP、原始p值和FDR校正值
data_VIP_P <- data.frame(
  VIP = data_VIP,
  P_value = p_values,
  FDR = fdr,
  effect_size = effect_sizes,
  log2FC = log2fc_values
)

# 将结果合并到原始数据
data_VIP_P_data <- cbind(raw_data, data_VIP_P)
data_VIP_P_data$name = rownames(data_VIP_P_data)
head(data_VIP_P_data)

# 保存完整结果
write.csv(data_VIP_P_data, "OPLSDA_VIP_P_JZL_CK_complete.csv", row.names = FALSE)

# 筛选显著特征：VIP > 1 且 FDR < 0.05（或0.01）
data_VIP_P_select <- data_VIP_P_data[data_VIP_P_data$VIP > 1 & data_VIP_P_data$P_value < 0.05 & abs(data_VIP_P_data$effect_size) > 0.8 & abs(data_VIP_P_data$log2FC) > 0.5,]
# 或者更严格的筛选：data_VIP_P_select <- data_VIP_P_data[data_VIP_P_data$VIP > 1 & data_VIP_P_data$FDR < 0.01, ]

# 查看筛选结果数量
cat("筛选到的显著特征数量:", nrow(data_VIP_P_select), "\n")
write.csv(data_VIP_P_select, "OPLSDA_VIP_P_JZL_CK_select.csv", row.names = FALSE)
