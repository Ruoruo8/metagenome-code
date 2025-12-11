rm(list=ls())#clear Global Environment
#设置工作目录
setwd("D:/metagenome/moshpit_tutorial/results/09α")
#安装包
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggsignif")
#install.packages("vegan")
#install.packages("ggprism")
#install.packages("picante")
#install.packages("dplyr")
#install.packages("RColorBrewer")
#加载包
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)
##导入数据，所需是数据行名为样本名、列名为OTUxxx的数据表
df <- read.csv("Genus_fraction.csv",header = T, row.names = 1, sep = ",")
#使用vegan包计算多样性指数
Shannon <- diversity(df, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(df, index = "simpson", MARGIN = 2, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)#spe.rich =sobs
###将以上多样性指数统计成表格
index <- as.data.frame(cbind(Shannon, Simpson, Richness))

#计算chao，ace指数
Chao1 <- function(x){
  S.obs <- sum(x  != 0)
  F1 <- sum(x  == 1)
  F2 <- sum(x  == 2)
  chao1 <- S.obs + F1*(F1-1)/(2*(F2+1))
  return(chao1)
}

ACE <- function(x, rare_threshold = 10){
  bord.num <- length(table(x))
  rare_threshold <- ifelse(rare_threshold < bord.num-1, rare_threshold, bord.num-1)
  s.rare <- sum(x <= rare_threshold & x > 0)
  singles <- sum(x == 1)
  s.abund <- sum(x > rare_threshold)
  freq_counts <- data.frame(table(x))[1:(rare_threshold+1),]
  nums <- seq(0,rare_threshold)
  N.rare.gamma <- sum(freq_counts[,2] * nums *  (nums-1))
  N.rare.nogamma <- sum(freq_counts[,2] * nums)
  c_ace = 1 - singles / N.rare.nogamma
  top <- s.rare * N.rare.gamma
  bottom <- c_ace * N.rare.nogamma * (N.rare.nogamma - 1)
  gamma_ace = (top / bottom) - 1  
  if (gamma_ace < 0){
    gamma_ace = 0
  }
  r <- s.abund + (s.rare / c_ace) + ((singles / c_ace) * gamma_ace)
  return(r)
}

# 初始化结果矩阵
result <- data.frame(
  Sample = colnames(df),
  Chao1 = numeric(ncol(df)),
  ACE = numeric(ncol(df))
)

# 计算每个样本的 Chao1 和 ACE
for (i in 1:ncol(df)) {
  sample_data <- df[, i]
  result$Chao1[i] <- Chao1(sample_data)
  result$ACE[i] <- ACE(sample_data)
}

#将obs，chao，ace指数与前面指数计算结果进行合并
index$Chao <- result[,2]
index$Ace <- result[,3]
#计算Pielou及覆盖度
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)
#导出表格
write.table(cbind(sample=c(rownames(index)),index),'diversity.index_fraction.txt', row.names = F, sep = '\t', quote = F)

#读入文件
index <- read.delim('diversity.index_fraction.txt', header = T, row.names = 1)
index$samples <- rownames(index)#将样本名写到文件中
#读入分组文件
groups <- read.csv('group.csv',header = T)
colnames(groups)[1:2] <- c('samples','group')#改列名
#合并分组信息与多样性指数
df2 <- merge(index,groups,by = 'samples')

#Shannon
p1 <- ggplot(df2,aes(x=group,y=Shannon))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+ #图例位置
  scale_fill_manual(values = c("#A5C496","#C7988C","#8891DB"))+  #指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45) # 可选值有 0，45，90，270
p1
ggsave("shannon_fraction.pdf",p1,width = 6, height = 4)

#Simpson
p2 <- ggplot(df2,aes(x=group,y=Simpson))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  scale_fill_manual(values = c("#A5C496","#C7988C","#8891DB"))+  #指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45) # 可选值有 0，45，90，270
p2
ggsave("Simpson_fraction.pdf",p2,width = 6, height = 4)

#Ace
p3 <- ggplot(df2,aes(x=group,y=Ace))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  scale_fill_manual(values = c("#A5C496","#C7988C","#8891DB"))+  #指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45) # 可选值有 0，45，90，270
p3
ggsave("Ace_fraction.pdf",p3,width = 6, height = 4)

#Chao
p4 <- ggplot(df2,aes(x=group,y=Chao))+#指定数据
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+#添加误差线,注意位置，放到最后则这条线不会被箱体覆盖
  geom_boxplot(aes(fill=group), #绘制箱线图函数
               outlier.colour="white",size=0.8)+#异常点去除
  theme(panel.background =element_blank(), #背景
        axis.line=element_line(),#坐标轴的线设为显示
        plot.title = element_text(size=14))+#图例位置
  scale_fill_manual(values = c("#A5C496","#C7988C","#8891DB"))+  #指定颜色
  geom_jitter(width = 0.2)+#添加抖动点
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.8, # 坐标轴的粗细
              axis_text_angle = 45) # 可选值有 0，45，90，270
p4
ggsave("Chao_fraction.pdf",p4,width = 6, height = 4)

library("gridExtra")
library("cowplot")
plot_grid(p1,p2,p3,p4, labels=c('A','B','C','D'), ncol=2, nrow=2)#拼图及标注
ggsave('α多样性_fraction.pdf',width=12,height = 13)
