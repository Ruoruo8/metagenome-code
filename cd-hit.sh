# 合并多个样本的基因序列（蛋白、DNA）
cat sample1_prot.faa sample2_prot.faa ... >  prot.faa
cat sample1_nucl.fna sample2_nucl.fna ... >  nucl.fna
# cd-hit 运行构建非冗余基因集（蛋白）
cd-hit -i all_prot.faa -o all_prot_nonerude.faa -c 0.95 -T 8 -n 5 -d 0 -aS 0.9 -g 1 -sc 1 -sf 1 -M 0
# 获取了非冗余基因名称列表，再通过seqtk把对应的非冗余基因集（DNA）序列挑选出来
grep '>' all_prot_nonerude.faa | awk -F ' ' '{print $1}'|sed 's/>//g' > all_prot_nonerude.list
seqtk subseq all_nucl.fna all_prot_nonerude.list > all_nucl_nonerude.fna
seqkit replace -p "(.*)" -r "Unigene.{nr}" all_nucl_nonerude.fna > nucl_nonerude.fna
seqkit replace -p "(.*)" -r "Unigene.{nr}" all_prot_nonerude.faa > prot_nonerude.faa
# 构建非冗余基因集的bwa索引
bwa index nucl_nonerude.fna -p geneset_bwa
bioawk -c fastx '{print $name, length($seq)}' nucl_nonerude.fna > geneset_length.txt

