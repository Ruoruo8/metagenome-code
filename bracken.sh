# 运行bracken
for file in *.kraken; do
  # 获取样本名（去掉路径和 "_clean_R1.fastq" 部分）
  sample_name=$(basename "$file" ".kraken")
  # 执行命令
  bracken -d ./k2db/ -i "${sample_name}.kreport" -o "${sample_name}.bracken.G" -l G -t 16
done

# 分类水平，如 D,P,C,O,F,G,S