# =========================
# 0. 加载 R 包
# =========================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ALDEx2")) BiocManager::install("ALDEx2")

library(ALDEx2)
library(dplyr)
library(tibble)
library(phyloseq)
library(ggplot2)
library(meta)

# =========================
# 1. 读取输入数据
# =========================
abund <- read.csv("Genus_counts.csv", header = TRUE, row.names = 1, check.names = FALSE)
group <- read.csv("group.csv", header = TRUE)

group$samples <- as.character(group$samples)
group$Group  <- as.factor(group$Group)

# 匹配样本顺序
abund <- abund[, group$samples]

# =========================
# 2. 自动生成 3 组的 pairwise 组合
# =========================
groups <- unique(as.character(group$Group))
comparisons <- combn(groups, 2, simplify = FALSE)   # 所有组合(A vs B, A vs C, B vs C)

# 输出文件列表保存
result_list <- list()

# =========================
# 3. 对每个 pairwise 运行 ALDEx2
# =========================
for (comp in comparisons) {
  
  g1 <- comp[1]
  g2 <- comp[2]
  message("Running ALDEx2: ", g1, " vs ", g2)
  
  # 取两组的样本
  idx <- which(group$Group %in% c(g1, g2))
  abund_sub  <- abund[, idx]
  group_sub  <- as.character(droplevels(group$Group[idx]))
  
  # ALDEx2 分析
  aldex_res <- aldex.clr(
    abund_sub,
    conds = group_sub,
    mc.samples = 128,
    denom = "all"
  )
  
  aldex_tt <- aldex.ttest(aldex_res)
  aldex_eff <- aldex.effect(aldex_res)
  
  # 合并结果
  aldex_df <- cbind(
    taxon = rownames(aldex_tt),
    aldex_tt,
    aldex_eff
  ) %>%
    as.data.frame() %>%
    mutate(
      Comparison = paste(g1, "vs", g2),
      BH_FDR = p.adjust(wi.ep, method = "BH")  # 默认使用 wilcox p-value 做 FDR
    )
  aldex_df$FDR_Qvalue <- p.adjust(aldex_df$we.ep, method = "BH")
  result_list[[paste(g1, "_vs_", g2, sep="")]] <- aldex_df
  
  # 保存结果文件
  write.csv(aldex_df,
            paste0("ALDEx2_", g1, "_vs_", g2, "_results_counts.csv"),
            row.names = FALSE)
}

# =========================
# 4. 合并所有结果
# =========================
all_results <- dplyr::bind_rows(result_list)
write.csv(all_results, "ALDEx2_all_pairwise_results.csv", row.names = FALSE)

message("ALDEx2 全部 pairwise 分析完成！")
