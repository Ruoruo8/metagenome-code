for file in *.contigs.fa; do
	file_name=$(basename "$file" .contigs.fa)
	seqkit seq -m 500 "$file" --remove-gaps > "${file_name}.contigs_500.fa"
	sed -i "s/>/>${file_name}_/g" "${file_name}.contigs_500.fa"
done
