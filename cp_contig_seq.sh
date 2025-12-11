#!/bin/bash

# 定义主文件夹路径
main_dir="/mnt/d/metagenome/03megahit"

# 遍历clean_data文件夹中的子文件夹
for subfolder in "$main_dir"/*; do
    # 确保只处理目录
    if [ -d "$subfolder" ]; then
        # 提取子文件夹名称
        subfolder_name=$(basename "$subfolder")

        # 遍历子文件夹中的所有文件
        for file in "$subfolder"/final.contigs.fa; do
            if [ -f "$file" ]; then
                # 构造目标文件路径
                target_file="${subfolder_name}.contigs.fa"
                # 复制文件并重命名
                cp "$file" "$target_file"
            fi
        done
    fi
done
