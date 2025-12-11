#!/bin/bash

# 当前目录为 metaphlan，定义输入数据文件夹
input_folder="../clean_data"

# 遍历 clean_data 文件夹中的所有 R1 文件
for file_r1 in "$input_folder"/*_clean_R1.fastq.gz; do
  # 获取样本名（去掉路径和 "_clean_R1.fastq" 部分）
  sample_name=$(basename "$file_r1" "_clean_R1.fastq.gz")

  # 构造 R2 文件路径
  file_r2="$input_folder/${sample_name}_clean_R2.fastq.gz"

  # 确保 R2 文件存在
  if [[ -f "$file_r2" ]]; then
    echo "Processing sample: $sample_name"

    # 构造输出文件名
    bowtie2_out="${sample_name}_bowtie2.bz2"
    metaphlan_out="${sample_name}_metaphlan.tsv"

    # 执行命令
    metaphlan "$file_r1","$file_r2" --input_type fastq --bowtie2out "$bowtie2_out" --output_file "$metaphlan_out" --nproc 8
  else
    echo "Warning: R2 file missing for sample: $sample_name"
  fi
done
