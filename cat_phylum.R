# 读取数据，不将第一列设为行名（因为会有重复）
asv_data <- read.csv("phylum.csv", header=TRUE, sep=",", check.names=FALSE)

# 查看数据结构
cat("数据维度:", dim(asv_data), "\n")
cat("列名:", colnames(asv_data), "\n")
head(asv_data)

# 方法1: 使用dplyr按Phylum合并计数
library(dplyr)

# 假设第一列名为"Phylum"（如果不是，请调整）
# 如果第一列没有列名，先设置一个
if (colnames(asv_data)[1] == "") {
  colnames(asv_data)[1] <- "Phylum"
}

phylum_summary <- asv_data %>%
  group_by(Phylum) %>%
  summarise(across(everything(), sum))

# 查看结果
cat("合并后的维度:", dim(phylum_summary), "\n")
head(phylum_summary)

# 保存结果
write.table(phylum_summary, "phylum_summary_table.txt", sep="\t", row.names=FALSE, quote=FALSE)