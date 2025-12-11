#!/bin/bash

total=0
count=0

# 遍历当前目录下所有.faa文件
for file in *.faa; do
	    # 跳过不匹配的情况
	        [ -e "$file" ] || continue
		    
		    # 统计单个文件的基因数量
		        genes=$(grep -c "^>" "$file")
			    echo "$file: $genes"
			        
			        # 累加总数
				    total=$((total + genes))
				        count=$((count + 1))
				done

				# 结果输出
				if [ $count -eq 0 ]; then
					    echo "未找到任何.faa文件"
				    else
					        echo "------------------------"
						    echo "总文件数：$count"
						        echo "总基因数：$total"
				fi
