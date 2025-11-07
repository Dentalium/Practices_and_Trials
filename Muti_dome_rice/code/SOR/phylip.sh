# 批量命令行执行phylip，动态配置配置文件

for no in {1..72}
do
	echo ${no}
	
	echo Fst_matrix/${no}.phylip >> phylip.config
	echo O >> phylip.config
	echo 5 >> phylip.config
	echo L >> phylip.config
	echo Y >> phylip.config

	~/biosoft/phylip-3.697/exe/neighbor < phylip.config > tree/${no}.log 2>&1

	mv outfile tree/${no}.neighbor.out
	mv outtree tree/${no}.nwk

	rm phylip.config
done
