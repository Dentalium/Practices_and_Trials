output_file="tree/combined_nwk.txt"

touch ${output_file}

for i in {1..72}
do
	echo -ne "${i}\t" >> ${output_file} 
	cat tree/${i}.nwk | tr -d '\n' >> ${output_file}
	echo "" >> ${output_file}
done
