# 读取PSR区间，计算每个区间在亚群之间的总FST
# vcftools不适合求任意区间平均FST，通过读取日志文件克服这个困难

PSR_list="../PSR_2.txt"
IFS=$'\n'

# 定义群体列表
pops=("indica" "japonica" "niv1" "niv2" "outer" "ruf1" "ruf2")

cat ${PSR_list} | while read psr
do
	no=$(echo ${psr} | cut -d $'\t' -f2)

	echo No ${no}
	
	mkdir -p Fst/${no}

	# 嵌套循环，遍历所有两两组合
	for i in ${!pops[@]}; do
		for j in ${!pops[@]}; do
			# 确保 i < j
			if [ $i -lt $j ]; then
				pop1=${pops[$i]}
				pop2=${pops[$j]}
				echo "$pop1 vs $pop2"
			       	
				# 运行 VCFtools
				vcftools --gzvcf PSR_extract/${no}.vcf.gz \
					--weir-fst-pop groups/${pop1}.txt \
					--weir-fst-pop groups/${pop2}.txt \
					--out Fst/${no}/${no}_${pop1}_${pop2}
			fi
		done
	done
done


