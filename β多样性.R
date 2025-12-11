#beta多样性
#数据载入
rm(list=ls())

#install.packages("ecodist")
library(vegan)
library(ecodist)
library(ggplot2)

otu <- read.csv("Genus_fraction.csv",header = T, row.names = 1, sep = ",")
totu <- t(otu)
spe <- as.data.frame(totu) 

#PCA
pca=rda(spe,scale=T) # scale=T，表示单位方差标准化，先中心化再均值化，物种数据均值为0，方差为1.
pca
## 绘图数据提取
s.pca=pca$CA$u # 提取样本特征值
s.pca
e.pca=pca$CA$v # 提取物种特征值
e.pca # 可执行选取排序轴绘制散点图
eig=pca$CA$eig # 提取特征根，计算PC1和PC2的解释度，横纵坐标标签

sample_scores <- as.data.frame(s.pca)
sample_scores$Sample <- rownames(sample_scores)

# 读取分组信息
sam <- read.csv('group.csv',row.names = 1)  # 第一列是样本名称，第二列是分组信息

# 将分组信息添加到样本坐标数据中
sample_scores$sample <- rownames(sample_scores)  # 添加样本名称列
sample_scores <- merge(sample_scores, sam, by.x = "sample", by.y = "row.names")  # 合并分组信息

# 提取特征根并计算解释度
explained_variance <- eig / sum(eig) * 100

pca_plot = ggplot(data = sample_scores,aes(x=PC1,y=PC2))+
  geom_point(aes(color=group,shape=group),size=3)+
  scale_colour_manual(values = c("#A5C496","#C7988C","#8891DB"))+
  stat_ellipse(aes(color = group), type = "t", level = 0.95) +  # 95% 置信区间
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")+
  # 添加样本名称
  geom_text(aes(label = sample), vjust = -0.5, size = 3) + 
  # 添加坐标轴标签
  xlab(paste("PC1 (", round(explained_variance[1], 2), "%)")) +
  ylab(paste("PC2 (", round(explained_variance[2], 2), "%)")) 

# 显示图形
print(pca_plot)
ggsave("PCA_plot_fraction.pdf", pca_plot, width = 8, height = 6, dpi = 300)


#PCoA
# 使用vegdist()函数求样品间距离矩阵
spe.dist.bray<-vegdist(spe,method = 'bray') # 距离矩阵计算有多种方法可以选择，"manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord".默认是bray
spe.dist.bray
# 将距离矩阵进行主坐标分析
pcoa = cmdscale(spe.dist.bray, k=3, eig=TRUE) # k设置坐标轴数量，推荐用3或者样本数-1; eig=TRUE返回特征根，默认不输出特征根结果
pcoa
poi = pcoa$points
eigval = pcoa$eig
pcoa_eig = (pcoa$eig)[]/sum(pcoa$eig)*100  #每个PCoA的权重值，这里对应下文ggplot中的lab
poi = as.data.frame(poi)

# 将分组信息添加到样本坐标数据中
poi$sample <- rownames(poi)  # 添加样本名称列
poi <- merge(poi, sam, by.x = "sample", by.y = "row.names")  # 合并分组信息

pcoa_plot = ggplot(data = poi,aes(x=V1,y=V2))+
  geom_point(aes(color=group,shape=group),size=3)+
  scale_colour_manual(values = c("#A5C496","#C7988C","#8891DB"))+
  stat_ellipse(aes(color = group), type = "t", level = 0.95) +  # 95% 置信区间
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")+
  # 添加样本名称
  geom_text(aes(label = sample), vjust = -0.5, size = 3) + 
  # 添加坐标轴标签
  xlab(paste("PCoA1 (", round(pcoa_eig[1], 2), "%)")) +
  ylab(paste("PCoA2 (", round(pcoa_eig[2], 2), "%)")) 
# 显示图形
print(pcoa_plot)
ggsave("PCoA_plot_fraction.pdf", pcoa_plot, width = 8, height = 6, dpi = 300)

# NMDS
#distance.bray<-vegdist(spe,method = 'bray')
#NMDS<-metaMDS(distance.bray,k=3)
NMDS<-metaMDS(spe,distance = "bray", k=3) # k设置维度数
NMDS
# 提取样本坐标数据
NMDS_poi = NMDS$points
# NMDS胁迫系数stress：检验NMDS分析结果的优劣，stress<0.2,可用NMDS的二维点图表示，其图形有一定解释意义；stress<0.1,认为是一个好的排序；stress<0.05,具有很好的代表性。
NMDS_stress = NMDS$stress
NMDS_poi = as.data.frame(NMDS_poi)
# 将分组信息添加到样本坐标数据中
NMDS_poi$sample <- rownames(NMDS_poi)  # 添加样本名称列
NMDS_poi <- merge(NMDS_poi, sam, by.x = "sample", by.y = "row.names")  # 合并分组信息

NMDS_plot = ggplot(data = NMDS_poi,aes(x=MDS1,y=MDS2))+
  geom_point(aes(color=group,shape=group),size=3)+
  scale_colour_manual(values = c("#A5C496","#C7988C","#8891DB"))+
  stat_ellipse(aes(color = group), type = "t", level = 0.95) +  # 95% 置信区间
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(aes(yintercept=0), colour="#BEBEBE", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="#BEBEBE", linetype="dashed")+
  # 添加样本名称
  geom_text(aes(label = sample), vjust = -0.5, size = 3) + 
  # 添加坐标轴标签
  xlab(paste("NMDS1")) +
  ylab(paste("NMDS2"))
# 显示图形
print(NMDS_plot)
ggsave("NMDS_plot_fraction.pdf", NMDS_plot, width = 8, height = 6, dpi = 300)