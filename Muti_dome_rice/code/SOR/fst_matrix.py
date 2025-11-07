def extract_wmFst(input_fst):
    '''
    从vcf日志中提取加权平均Fst，转换为小数点后5位
    '''
    with open(input_fst) as f:
        for i in f:
            if i.startswith("Weir and Cockerham weighted Fst estimate"):
                wmFst=float(i.rstrip().split(sep=":")[1].lstrip(" "))
                wmFst=f"{wmFst:.5f}"
                break
    
    return wmFst




pop_list=["indica","japonica","niv1","niv2","outer","ruf1","ruf2"]

for no in range(1, 73):
    output_matrix=f"Fst_matrix/{no}.phylip"

    with open(output_matrix, "w") as f:
        # header
        f.write(" 7\n")

        for i in range(7):
            fst_list=[]
            pop1=pop_list[i]
            # 样本名
            f.write(f"{pop1:<10}")
            for j in range(i):
                pop2=pop_list[j]

                input_fst=f"Fst/{no}/{no}_{pop2}_{pop1}.log"
                fst_list.append(extract_wmFst(input_fst=input_fst))
            
            f.write(" ".join(fst_list))
            f.write("\n")

