#!/usr/bin/env python3

"""
模仿ZhangFM的perl脚本，滑窗统计群体遗传学统计量，有随机访问能力，提高效率
略微修改了原脚本分组名称
原脚本header多了segment一列，但是在脚本中不会输出这一列，因此将其删除
作者计算theta和tajima's D时定义了“等效分离位点数”用于处理复等位基因位点，这个是否正确？
"""
import sys
import math

# 定义计算群体遗传统计量的函数
# 所有函数均传入对应样本的一个窗口内的所有位点的field，列表
def cal_k(fr):
    """计算未归一化的pi，即区间杂合度的总和"""
    k = 0.0
    # 遍历所有位点
    for freq_str in fr:
        frqs = [float(x) for x in freq_str.split(':')] # 各等位基因频数
        AN = len(frqs) # 等位基因数量，最小为2
        mcross = 0.0
        # 所有等位基因对的频数组合乘积之和（交叉项）
        for i in range(AN-1):
            for j in range(i+1, AN):
                mcross += frqs[i] * frqs[j]
        
        AC = sum(frqs) # 等位基因总数
        k += 2 * mcross / (AC * (AC - 1))

    return k

def cal_dxy(p1, p2, ln):
    """计算群体间平均核苷酸差异
    dxy=$\sum_{i\neq j}$
    p1、p2：对应样本的所有位点的field
    ln：区间长度
    """
    dxy = 0.0

    for i in range(len(p1)):
        frqs1 = [float(x) for x in p1[i].split(':')]
        frqs2 = [float(x) for x in p2[i].split(':')]
        
        mcross = 0.0
        AN = len(frqs1)
        AC1 = sum(frqs1)
        AC2 = sum(frqs2)
        
        # 计算不同等位基因的频数乘积之和
        for j in range(AN):
            for k in range(AN):
                if j == k:
                    continue
                else:
                    mcross += frqs1[j] * frqs2[k]

        dxy += mcross / (AC1 * AC2)
    
    dxy = dxy / ln # 按区间长度平均
    return dxy

def cal_ss(fr):
    """区间内等效分离位点数量的总和
    每个位点的等效位点数等于等位基因数减1
    """
    ss = 0
    for freq_str in fr:
        frqs = [float(x) for x in freq_str.split(':')]
        # 统计频数不为0的等位基因数量
        i = sum(1 for x in frqs if x > 0)
        ss += i - 1
    return ss

def cal_sn(fr):
    """计算区间内平均AC，取整数"""
    n = 0
    fn = 0.0
    for freq_str in fr:
        frqs = [float(x) for x in freq_str.split(':')]
        i = sum(1 for x in frqs if x > 0)
        AC = sum(frqs)  # 总等位基因计数
        if i - 1 > 0:   # 只统计分离位点
            n += 1
            fn += AC
    sn = 0
    if n > 0:
        sn = int(fn / n + 0.5)
    return sn

def cal_theta(ss, sn, ln):
    """计算Watterson's θ
    """
    theta = 0.0

    a = 0.0
    for i in range(1, sn):
        a += 1.0 / i
    theta = ss / a / ln

    return theta

def cal_tajimaD(k, ss, sn):
    """计算Tajima's D"""
    td = 0.0
    if ss > 1 and sn > 0:
        a1 = 0.0
        for i in range(1, sn):
            a1 += 1.0 / i
        
        b1 = (sn + 1) / (3.0 * (sn - 1))
        c1 = b1 - 1.0 / a1
        e1 = c1 / a1
        
        a2 = 0.0
        for i in range(1, sn):
            a2 += 1.0 / (i * i)
        
        b2 = 2.0 * (sn * sn + sn + 3) / (9.0 * sn * (sn - 1))
        c2 = b2 - (sn + 2) / (a1 * sn) + a2 / (a1 * a1)
        e2 = c2 / (a1 * a1 + a2)
        
        td = (k - ss / a1) / math.sqrt(e1 * ss + e2 * ss * (ss - 1))
    
    return td

def main():
    if len(sys.argv) != 6:
        print("Usage: python script.py <frq> <frq_idx> <wd_size> <wd_step> <out>")
        sys.exit(1)
    
    snp_frq = sys.argv[1]   # 分窗统计的等位基因频数
    frq_idx = sys.argv[2]   # 等位基因频数文件的索引
    wd_size = int(sys.argv[3]) * 1000  # 滑窗大小 (kp)
    wd_step = int(sys.argv[4]) * 1000  # 滑窗步长 (kp)
    out_file = sys.argv[5]  # output file
    
    chromosomes = [str(i) for i in range(1,13)]
    pop_names = ['dome', 'wild']
    
    # 读取索引，窗口代号:指针位置
    pos_idx = {}
    with open(frq_idx) as f:
        for line in f:
            fields = line.rstrip().split('\t')
            pos_idx[fields[0]] = int(fields[1])
    
    # 打开输出文件
    with open(out_file, 'w') as out_f:
        # 表头
        header = ['chr', 'start_pos', 'end_pos', 'snp_num']
        for pop in pop_names:
            header.extend([f"{pop}-pi(1kb)", f"{pop}-theta-W(1kb)", f"{pop}-Tajima-D/sample"])
        header.extend([f"{pop_names[0]}-{pop_names[1]}-dxy(1kb)", 
                    f"{pop_names[0]}-{pop_names[1]}-fst", 
                    f"{pop_names[0]}-{pop_names[1]}-rod"])
        out_f.write('\t'.join(header) + '\n')
        
        # 打开等位基因频数文件
        with open(snp_frq, 'r') as in_f:
            # 读取（跳过）表头
            header_line = in_f.readline()
            
            # 分组信息的索引（下标），用来在统计时区分两个分组
            pop_pos = {"dome":5, "wild":6}
            
            # 滑窗计算统计量
            # 遍历染色体
            for chromosome in chromosomes:
                # 滑窗起点
                start_pos = 1
                flag = True
                
                print("chr:",chromosome)

                while flag:
                    end_pos = start_pos + wd_size - 1
                    length = end_pos - start_pos + 1 # 滑窗实际大小（恒等于wd_size，意义不明）
                    
                    print(start_pos)

                    # 准备输出行
                    output_fields = [chromosome, str(start_pos), str(end_pos)]
                    
                    # 初始化群体频率字典
                    # key为分组名称，值为列表，储存对应当前窗口所有位点的fields
                    pop_frq = {pop: [] for pop in pop_names}
                    
                    # 滑窗序号
                    # 注意是索引中定义的滑窗而不是统计的滑窗，统计的滑窗是100kb（10倍）
                    pos_marker = start_pos // 10000
                    # 这里模仿了原perl脚本中的“找存在的最小窗口序号”的逻辑，以防索引不完整，虽然似乎没什么用
                    # 如果一直找不到对应窗口，索引变成-1（或者死循环），后面移动指针时就会报错
                    while f"{chromosome}-{pos_marker}" not in pos_idx and pos_marker >= 0:
                        pos_marker -= 1
                    
                    # 指针移动到指定地点，准备读取
                    if pos_marker >= 0:
                        in_f.seek(pos_idx[f"{chromosome}-{pos_marker}"])
                    else:
                        print("Index not integrate!")
                        exit(1)
                    
                    # 读取窗口内的SNP
                    snp_count = 0
                    for line in in_f:
                        fields = line.rstrip().split(sep='\t')
                        
                        chr_name = fields[0]
                        pos = int(fields[1])
                        
                        # 检查是否在当前窗口内
                        if chr_name == chromosome and start_pos <= pos <= end_pos:
                            snp_count += 1
                            for pop in pop_names:
                                pop_frq[pop].append(fields[pop_pos[pop]])
                        
                        # 如果超过当前窗口或染色体，停止读取
                        elif chr_name != chromosome or pos > end_pos:
                            # 如果是不同染色体或文件结束，设置flag=False，换下一条染色体
                            if chr_name != chromosome or not line:
                                flag = False
                            break
                    
                    # 输出SNP计数
                    output_fields.append(str(snp_count))
                    
                    # 计算统计量
                    k = {}; ss = {}; sn = {}; pi = {}; th = {}; tj = {}
                    
                    for pop in pop_names:
                        k[pop] = cal_k(pop_frq[pop])
                        pi[pop] = 1000 * k[pop] / length if length > 0 else 0

                        ss[pop] = cal_ss(pop_frq[pop])
                        sn[pop] = cal_sn(pop_frq[pop])
                        th[pop] = 1000 * cal_theta(ss[pop], sn[pop], length) if length > 0 else 0
                        tj[pop] = cal_tajimaD(k[pop], ss[pop], sn[pop])
                        
                        output_fields.extend([f"{pi[pop]:10.5f}", 
                                            f"{th[pop]:10.5f}", 
                                            f"{tj[pop]:10.5f}/{sn[pop]}"])
                    
                    # 计算群体间统计量
                    p1, p2 = pop_names[0], pop_names[1]
                    dxy = cal_dxy(pop_frq[p1], pop_frq[p2], length) if length > 0 else 0
                    dxy = 1000 * dxy
                    
                    if dxy == 0:
                        fst = 0
                    else:
                        fst = (2 * dxy - pi[p1] - pi[p2]) / (2 * dxy)
                        if fst < 0:
                            fst = 0
                    
                    if pi[p1] == 0:
                        rod = 999 # 原脚本如此
                    else:
                        rod = pi[p2] / pi[p1]
                    
                    output_fields.extend([f"{dxy:10.5f}", f"{fst:10.5f}", f"{rod:10.5f}"])
                    
                    # 写入输出
                    out_f.write('\t'.join(output_fields) + '\n')
                    
                    # 移动到下一个窗口
                    start_pos += wd_step



if __name__ == "__main__":
    main()