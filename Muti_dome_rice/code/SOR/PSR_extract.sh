PSR_list="../PSR_2.txt"
vcf_input="../../../1355_wild+dome_maj.vcf.gz"

IFS=$'\n'

cat ${PSR_list} | while read i
do
	chr=$(echo ${i} | cut -d $'\t' -f1)
	no=$(echo ${i} | cut -d $'\t' -f2)
	start=$(echo ${i} | cut -d $'\t' -f3)
	end=$(echo ${i} | cut -d $'\t' -f4)
	
	echo ${no}

	bcftools view -r ${chr}:${start}-${end} -o PSR_extract/${no}.vcf.gz ${vcf_input}
done
