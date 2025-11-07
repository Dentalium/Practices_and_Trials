'''
这个脚本用来分组统计等位基因频数
'''

import pysam
import argparse
import os

parser=argparse.ArgumentParser()

parser.add_argument("-S", required=True, help="The grouped list of samples")
parser.add_argument("-v", required=True, help="Input vcf.gz file")

args=parser.parse_args()

sample_list=args.S
input_vcf=args.v

output_file=os.path.basename(input_vcf.rstrip("vcf.gz"))+"freq.txt"

def freq_stat(this_sample:pysam.VariantRecordSample, n_alleles:int)->tuple:
    '''
    统计一个样本一个位点上各等位基因频数
    返回一个元组，元素为等位基因0，1，2，...的数量
    '''
    this_GT=this_sample["GT"]

    if None in this_GT:
        return tuple(0 for i in range(n_alleles))
    else:
        return tuple(this_GT.count(i) for i in range(n_alleles))

# 解析样品列表
samples_wild=[]
samples_dome=[]

with open(sample_list) as f:
    for i in f:
        this_line=i.rstrip().split(sep="\t")
        this_sample=this_line[0]
        this_sample_type=this_line[1]
        if this_sample_type=="wild":
            samples_wild.append(this_sample)
        elif this_sample_type=="dome":
            samples_dome.append(this_sample)
        else:
            print("Sample list error!")
            exit(1)

with open(output_file, "w") as f:
    f.write("#CHROM\tPOS\tREF\tALT\tall\tdome\twild\n")

    # 读取vcf文件
    myvcf=pysam.VariantFile(input_vcf, mode="r")
    try:
        for rec in myvcf.fetch():
            # 含有*的位点跳过
            if "*" in rec.alts:
                continue

            alts_str=",".join(list(rec.alts))
            n_alleles=len(rec.alleles)

            # 统计wild样品
            freq_wild=[0 for j in range(n_alleles)]
            for i in samples_wild:
                freq_this_sample=freq_stat(this_sample=rec.samples[i], n_alleles=n_alleles)
                # 累加
                for j in range(n_alleles):
                    freq_wild[j]+=freq_this_sample[j]
            
            # 统计dome样品
            freq_dome=[0 for j in range(n_alleles)]
            for i in samples_dome:
                freq_this_sample=freq_stat(this_sample=rec.samples[i], n_alleles=n_alleles)
                # 累加
                for j in range(n_alleles):
                    freq_dome[j]+=freq_this_sample[j]
            
            # 求所有样品之和
            freq_all=[]
            for i in range(n_alleles):
                freq_all.append(freq_wild[i]+freq_dome[i])
            
            # 转换格式
            freq_all_str=":".join([str(i) for i in freq_all])
            freq_dome_str=":".join([str(i) for i in freq_dome])
            freq_wild_str=":".join([str(i) for i in freq_wild])
            
            # 输出
            f.write(f"{rec.chrom}\t{rec.pos}\t{rec.ref}\t{alts_str}\t{freq_all_str}\t{freq_dome_str}\t{freq_wild_str}\n")

    finally:
        myvcf.close()
