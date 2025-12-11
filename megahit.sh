#!/bin/bash

# 当前目录为 metaphlan，定义输入数据文件夹
input_folder="../01clean_data"
# 定义临时文件夹，设置为 /home/ruoruo（或其他合适的路径）
temp_folder="/home/ruoruo/temp_megahit"

# 遍历 clean_data 文件夹中的所有 R1 文件
for file_r1 in "$input_folder"/*_clean_R1.fastq.gz; do
  # 获取样本名（去掉路径和 "_clean_R1.fastq" 部分）
  sample_name=$(basename "$file_r1" "_clean_R1.fastq.gz")

  # 构造 R2 文件路径
  file_r2="$input_folder/${sample_name}_clean_R2.fastq.gz"

  # 确保 R2 文件存在
  if [[ -f "$file_r2" ]]; then
    echo "Processing sample: $sample_name"

    # 执行 megahit 命令，指定临时文件夹作为输出目录
    megahit -1 "$file_r1" -2 "$file_r2" -o "$temp_folder/$sample_name"

    # 将结果文件移动到当前目录
    mv "$temp_folder/$sample_name" ./

    # 删除临时文件夹，清理中间文件
    rm -rf "$temp_folder/$sample_name"

    echo "Sample $sample_name processing completed and cleaned up."
  else
    echo "Warning: R2 file missing for sample: $sample_name"
  fi
done
