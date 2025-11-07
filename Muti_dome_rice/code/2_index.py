import sys
import os

'''
为分组统计的等位基因频数结果文件创建字节索引，以实现随机访问
'''

# 窗口宽度
sw_wide=10000

def create_genome_index(input_file, index_file):
    # 索引
    indices = {}
    # 窗口编号
    sw_keys = []

    current_chr = "chr0"
    pos_value = 0
    
    with open(input_file, 'rb') as f:  # 二进制模式读取
        line = f.readline()

        # 跳过header
        while line and line.startswith(b'#'):
            pos_value = f.tell()
            line = f.readline()
        
        # 处理数据行
        while line:
            fields = line.decode().rstrip().split('\t')

            chrom = fields[0]
            pos = int(fields[1])
            
            # 检查染色体变化
            if chrom != current_chr:
                current_chr = chrom
                #current_sw=-1
                current_sw=0

                # 根据原作者代码逻辑，必须在索引中写入染色体起点，否则一旦第一个窗口不存在SNP则会陷入死循环
                key = f"{chrom}-{0}"
                indices[key] = pos_value
                sw_keys.append(key)

            # 窗口编号
            sw_no=pos // sw_wide
            if sw_no != current_sw:
                current_sw=sw_no
                
                key = f"{chrom}-{sw_no}"
                indices[key] = pos_value
                sw_keys.append(key)
            
            pos_value = f.tell()
            line = f.readline()
    
    # 写入索引文件
    with open(index_file, 'w') as f:
        for key in sw_keys:
            f.write(f"{key}\t{indices[key]}\n")

if __name__ == "__main__":
    if sys.argv[1]=="-h":
        print(f"python {sys.argv[0]} <输入文件>")
        sys.exit(0)
    
    input_file = sys.argv[1]
    index_file = os.path.basename(input_file) + '.idx'
    
    # 创建索引（如果不存在）
    
    create_genome_index(input_file, index_file)