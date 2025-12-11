#!/bin/bash

# 定义输入数据文件夹
input_folder="../01clean_data"

# 遍历 clean_data 文件夹中的所有 R1 文件
for file_r1 in "$input_folder"/*_clean_R1.fastq.gz; do
  # 获取样本名（去掉路径和 "_clean_R1.fastq" 部分）
  sample_name=$(basename "$file_r1" "_clean_R1.fastq.gz")

  # 构造 R2 文件路径
  file_r2="$input_folder/${sample_name}_clean_R2.fastq.gz"

  # 确保 R2 文件存在
  if [[ -f "$file_r2" ]]; then
    echo "Processing sample: $sample_name"
    # 执行命令
    bwa mem -t 4 geneset_bwa "$file_r1" "$file_r2" | samtools view -bS - | samtools sort - > "${sample_name}_mapping_geneset.bam"
    samtools view -F 4 -F 256 -F 2048 "${sample_name}_mapping_geneset.bam" | awk '{if($3!="*") print $3}'|sort| uniq -c|awk 'BEGIN {FS=" ";OFS=","} {print $2,$1}' | awk 'BEGIN {FS=",";OFS=","} {if ($2 > 1) print $1"\\t"$2; else print $1"\\t0"}'|sed '1i gene\\tsample' > "${sample_name}.count"
  else
    echo "Warning: R2 file missing for sample: $sample_name"
  fi
done
