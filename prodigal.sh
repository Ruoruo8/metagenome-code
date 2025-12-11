#!/bin/bash

# 定义输入数据文件夹
input_folder="../04contig"

# 遍历 04contig 文件夹中的所有 .contigs_500.fa 文件
for file in "$input_folder"/*.contigs_500.fa; do
  # 获取样本名
  sample_name=$(basename "$file" ".contigs_500.fa")
  # 运行脚本
  prodigal -p meta -a "${sample_name}_prot.faa" -m -d "${sample_name}_nucl.fna" -o "${sample_name}_genes.gff" -f gff -s "${sample_name}.stat" -i "$file"
  
done
