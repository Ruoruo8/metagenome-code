import pandas as pd
import os

# 定义文件路径
file_dir = "/mnt/d/metagenome/moshpit_tutorial/results/06Unigenes/counts/"  # 替换为你的文件目录
output_file = "merged_counts.txt"  # 合并后的输出文件

# 初始化一个空的 DataFrame
merged_df = None

# 遍历目录中的所有 .count.txt 文件
for file_name in os.listdir(file_dir):
    if file_name.endswith(".count.txt"):
        file_path = os.path.join(file_dir, file_name)
        # 读取文件
        df = pd.read_csv(file_path, sep="\t")
        # 设置 gene 列为索引
        df.set_index("gene", inplace=True)
        # 重命名 sample 列为文件名
        df.rename(columns={"sample": file_name.replace(".count.txt", "")}, inplace=True)
        # 合并到主 DataFrame
        if merged_df is None:
            merged_df = df
        else:
            merged_df = merged_df.join(df, how="outer")

# 保存合并后的文件
merged_df.to_csv(output_file, sep="\t")
print(f"合并完成，结果已保存到 {output_file}")