# 安装 VennDiagram 包，如果已经安装过可以省略这一步
install.packages("VennDiagram")

# 加载 VennDiagram 包以使用其函数
library(VennDiagram)
# 导入数据
data <- read.csv("diff_venn.csv",header = T, sep=',')

venn.plot <- venn.diagram(
  x = list(Mcor_CK = data$Mcor_CK, Mmic_CK = data$Mmic_CK, Mmic_Mcor = data$Mmic_Mcor),
  filename = NULL,
  height = 600,
  width = 600,
  resolution = 300,
  col = "transparent",
  fill = c("#A5C496", "#C7988C", "#8891DB"),
  alpha = 0.5,
  cex = 0.45,
  cat.cex = 0.45
)
venn.plot
pdf("venn.pdf")
grid.draw(venn.plot)
dev.off()
