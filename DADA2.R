# 检查是否存在Biocondoctor安装工具，没有则安装
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",repo=site)
# 加载安装工具与安装DADA2，如果提示R包版本与R版本不匹配，可以根据提示信息选择安装合适的DADA2版本
library(BiocManager)
BiocManager::install("dada2")
# 加载DADA2包
library(dada2)
packageVersion("dada2")
path <- "D:/metagenome/16S/DADA2"
list.files(path)
# 返回测序正向文件完整文件名
fnFs <- sort(list.files(path, pattern= "_R1.fq", full.names = TRUE))
# 返回测序反向文件完整文件名
fnRs <- sort(list.files(path, pattern="_R2.fq", full.names = TRUE))
# 提取文件名中第一个`_`分隔的前文本作为样品名
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 检查提取出的文件名
sample.names
# 设置过滤文件的输出路径，将过滤后的文件存于\filtered
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq"))
# 过滤文件输出，统计结果保存于out
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=T)
head(out)
# 分别计算正向和反向序列错误率
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errR, nominalQ=TRUE)

# 去除正向序列数据中的重复序列
derepFs <- derepFastq(filtFs, verbose=TRUE)

# 去除反向序列数据中的重复序列
derepRs <- derepFastq(filtRs, verbose=TRUE)

# 基于错误模型进一步质控
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# 合并双端序列
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# 构建ASV表，amplicon sequence variant（ASV）表类似于我们传统的OTU表
seqtab <- makeSequenceTable(mergers)
# dim()的第二个值为扩增子序列个数，table(nchar())的统计结果表示每个读长下有多少个扩增子序列。
dim(seqtab)
# Inspect distribution of sequence lengths 查看序列长度分布
table(nchar(getSequences(seqtab)))

#去除嵌合体
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# 导出生成的ASV表（seqtab.nochim）
write.csv(seqtab.nochim,file="ASV.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

# 用训练器将序列分类，生成的文件将保存至工作目录，如果放在其它位置，请修改下面代码中path变量
taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "/silva_nr99_v138.2_toSpecies_trainset.fa.gz "), multithread=TRUE)
# 完成分类后，用参考序列数据包，对应填充数据信息
taxa <- addSpecies(taxa, paste0(path, "/silva_v138.2_assignSpecies.fa.gz"))
# 输出taxa文件
write.csv(taxa,file="taxa.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
# 另存物种注释变量，去除序列名，只显示物种信息
# Removing sequence rownames for display only
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
# 输出taxa.print文件
write.csv(taxa.print,file="taxa.print.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

# 安装DECIPHER，加载DECIPHER
BiocManager::install("DECIPHER")
library(DECIPHER)
packageVersion("DECIPHER")
# 转换ASV表为DNAString格式
# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
# 相关数据下载详见DECIPHER教程，并修改为下载目录
# CHANGE TO THE PATH OF YOUR TRAINING SET
load(paste0(path, "/SILVA_SSU_r138_2_2024.RData"))
# use all processors
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
# ranks of interest
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
# Removing sequence rownames for display only
taxa.print.DECIPHER <- taxa <- taxid
rownames(taxa.print.DECIPHER) <- NULL
head(taxa.print.DECIPHER)

# 导出注释结果文件taxa.print.DECIPHER
write.csv(taxa.print.DECIPHER,file="taxa.print.DECIPHER.CSV",append = FALSE, quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
