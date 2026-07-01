# Calculate the consistency in real cell lines data: HCT-116.

import os
import openpyxl

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.4 b/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

def GFHunter_RF(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'GFHunter_PacBio.csv')
    line = PacBio_file.readline()
    while line:
        if 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'GFHunter_cDNA_e50.csv')
    line = cDNA_file.readline()
    while line:
        if 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'GFHunter_dRNA.csv')
    line = dRNA_file.readline()
    while line:
        if 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def GFHunter_SF(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'GFHunter_PacBio.csv')
    line = PacBio_file.readline()
    while line:
        if 'flag = 5' in line or 'flag = 4' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'GFHunter_cDNA_e50.csv')
    line = cDNA_file.readline()
    while line:
        if 'flag = 5' in line or 'flag = 4' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'GFHunter_dRNA.csv')
    line = dRNA_file.readline()
    while line:
        if 'flag = 5' in line or 'flag = 4' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def GFHunter(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'GFHunter_PacBio.csv')
    line = PacBio_file.readline()
    while line:
        if 'flag' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'GFHunter_cDNA_e50.csv')
    line = cDNA_file.readline()
    while line:
        if 'flag' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'GFHunter_dRNA.csv')
    line = dRNA_file.readline()
    while line:
        if 'flag' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def Genion(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'Genion_PacBio.tsv')
    line = PacBio_file.readline()
    while line:
        gene1 = line.split('\t')[1].split('::')[0]
        gene2 = line.split('\t')[1].split('::')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'Genion_cDNA.tsv')
    line = cDNA_file.readline()
    while line:
        gene1 = line.split('\t')[1].split('::')[0]
        gene2 = line.split('\t')[1].split('::')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'Genion_dRNA.tsv')
    line = dRNA_file.readline()
    while line:
        gene1 = line.split('\t')[1].split('::')[0]
        gene2 = line.split('\t')[1].split('::')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def JAFFAL_HC(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'JAFFAL_PacBio.csv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        if 'HighConfidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'JAFFAL_cDNA.csv')
    line = cDNA_file.readline()
    line = cDNA_file.readline()
    while line:
        if 'HighConfidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'JAFFAL_dRNA.csv')
    line = dRNA_file.readline()
    line = dRNA_file.readline()
    while line:
        if 'HighConfidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def JAFFAL_C(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'JAFFAL_PacBio.csv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'JAFFAL_cDNA.csv')
    line = cDNA_file.readline()
    line = cDNA_file.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'JAFFAL_dRNA.csv')
    line = dRNA_file.readline()
    line = dRNA_file.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def JAFFAL(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'JAFFAL_PacBio.csv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        gene1 = line.split(',')[1].split(':')[0]
        gene2 = line.split(',')[1].split(':')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'JAFFAL_cDNA.csv')
    line = cDNA_file.readline()
    line = cDNA_file.readline()
    while line:
        gene1 = line.split(',')[1].split(':')[0]
        gene2 = line.split(',')[1].split(':')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'JAFFAL_dRNA.csv')
    line = dRNA_file.readline()
    line = dRNA_file.readline()
    while line:
        gene1 = line.split(',')[1].split(':')[0]
        gene2 = line.split(',')[1].split(':')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def LongGF(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'LongGF_PacBio.log')
    line = PacBio_file.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'LongGF_cDNA.log')
    line = cDNA_file.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'LongGF_dRNA.log')
    line = dRNA_file.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def fusionseeker(result_dir):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'fusionseeker_PacBio.txt')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        gene1 = line.split('\t')[1]
        gene2 = line.split('\t')[2]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'fusionseeker_cDNA.txt')
    line = cDNA_file.readline()
    line = cDNA_file.readline()
    while line:
        gene1 = line.split('\t')[1]
        gene2 = line.split('\t')[2]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l2.append(key)
        line = cDNA_file.readline()
    dRNA_file = open(result_dir + 'fusionseeker_dRNA.txt')
    line = dRNA_file.readline()
    line = dRNA_file.readline()
    while line:
        gene1 = line.split('\t')[1]
        gene2 = line.split('\t')[2]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l3.append(key)
        line = dRNA_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

# 新增：读取 ctat-LR-fusion 结果的函数 (支持 PacBio, cDNA, dRNA)
def ctat(result_dir):
    l1 = [] # PacBio
    l2 = [] # cDNA
    l3 = [] # dRNA
    
    # 1. 读取 PacBio 结果
    try:
        f = None
        for ext in ['_PacBio.tsv', '.PacBio.tsv', '_PB.tsv', '.PB.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                break
        if f:
            line = f.readline() # 跳过表头
            line = f.readline()
            while line:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    num = int(parts[1])
                    if num >= 2:
                        gene1 = parts[2]
                        gene2 = parts[5]
                        ls = [gene1, gene2]
                        ls.sort()
                        key = ls[0] + ':' + ls[1]
                        l1.append(key)
                line = f.readline()
            f.close()
    except FileNotFoundError:
        pass

    # 2. 读取 cDNA 结果
    try:
        f = None
        for ext in ['_cDNA.tsv', '.cDNA.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                break
        if f:
            line = f.readline()
            line = f.readline()
            while line:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    gene1 = parts[2]
                    gene2 = parts[5]
                    ls = [gene1, gene2]
                    ls.sort()
                    key = ls[0] + ':' + ls[1]
                    l2.append(key)
                line = f.readline()
            f.close()
    except FileNotFoundError:
        pass

    # 3. 读取 dRNA 结果
    try:
        f = None
        for ext in ['_dRNA.tsv', '.dRNA.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                break
        if f:
            line = f.readline()
            line = f.readline()
            while line:
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    num = int(parts[1])
                    if num >= 2:
                        gene1 = parts[2]
                        gene2 = parts[5]
                        ls = [gene1, gene2]
                        ls.sort()
                        key = ls[0] + ':' + ls[1]
                        l3.append(key)
                line = f.readline()
            f.close()
    except FileNotFoundError:
        pass

    new_l1 = list(set(l1))
    l1 = new_l1
    new_l2 = list(set(l2))
    l2 = new_l2
    new_l3 = list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    
    data = (l1, l2, l3)
    print(len(l1), len(l2), len(l3))
    return data

def JAFFA():
    l1 = []
    PacBio_file = open(path + '/result/NGS/JAFFA_MCF7.csv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        gene1 = line.split(',')[1].split(':')[0]
        gene2 = line.split(',')[1].split(':')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    l1.sort()
    
    return l1

def starfusion():
    l1 = []
    PacBio_file = open(path + '/result/NGS/starfusion_MCF7.tsv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        gene1 = line.split('\t')[0].split('--')[0]
        gene2 = line.split('\t')[0].split('--')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
    new_l1=list(set(l1))
    l1 = new_l1
    l1.sort()
    
    return l1

def compare(ngs, ls, r):
    result1 = []
    result2 = []
    for key in ngs:
        if key in ls:
            if key in r:
                result1.append(1)
                result2.append(0)
            else:
                result1.append(0)
                result2.append(1)
        else:
            result1.append(0)
            result2.append(0)
    return (result1, result2)

def read_result(dir):
    results = []
    file = open(dir + 'real_result.txt')
    line = file.readline()
    while line:
        g1 = line.split(':')[0]
        g2 = line.split('\n')[0].split(':')[1]
        name = line.split('\n')[0]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        results.append(key)
        line = file.readline()
    results = list(set(results))
    return results

if __name__ == "__main__":
    dir = path + '/result/MCF-7/'
    result = read_result(dir)
    gf_rf = GFHunter_RF(dir)
    gf_sf = GFHunter_SF(dir)
    gf = GFHunter(dir)
    lg = LongGF(dir)
    ja_hc = JAFFAL_HC(dir)
    ja_c = JAFFAL_C(dir)
    ja = JAFFAL(dir)
    fs = fusionseeker(dir)
    ge = Genion(dir)
    
    # 获取 ctat-LR-fusion 数据
    ct = ctat(dir)
    
    # 加入 ctat 的结果元组
    r = (gf_rf, gf_sf, gf, lg, ja_hc, ja_c, ja, fs, ge, ct)
    
    jan = JAFFA()
    sf = starfusion()
    ngs = []
    for key in jan:
        if key in sf:
            ngs.append(key)

    # ================= 新增：按照指定顺序排序 =================
    # 您提供的目标排序列表 (注意已转换为代码内部使用的 ':' 格式)
    custom_order = [
        "AHCYL1:RAD51C",
        "AK7:PAPOLA",
        "ARFGEF2:SULF2",
        "ATP1A1:ZFP64",
        "ATXN7L3:FAM171A2",
        "BCAS3:BCAS4",
        "BCAS4:ZMYND8",
        "CARM1:SMARCA4",
        "CCDC170:ESR1",
        "DEPDC1B:ELOVL7",
        "FCHO1:MYO9B",
        "GATAD2B:NUP210L",
        "GPR37L1:NAV1",
        "MATN2:POP1",
        "MYO6:SENP6",
        "NBPF6:SLC25A24",
        "PICALM:SYTL2",
        "RPS6KB1:VMP1",
        "SYAP1:TXLNG",
        "CA4:TANC2",
        "EIF4E2:GIGYF2",
        "EMCN:SMARCC1",
        "GPM6A:SEZ6L2",
        "TANC2:TLK2"
    ]
    
    # 使用自定义顺序对 ngs 列表排序
    # 如果 ngs 包含列表中有的项，按列表索引排序；如果有不在列表中的项，放到最后(9999)
    ngs.sort(key=lambda x: custom_order.index(x) if x in custom_order else 9999)
    # ==========================================================
            
    # 追加 ctat 的标题名
    ls = ['GFHunter RF', 'GFHunter RF + SF', 'GFHunter', 'LongGF', 'JAFFAL_HC', 'JAFFAL_C', 'JAFFAL', 'fusionseeker', 'Genion', 'ctat']
    list_types = ['PB', 'cDNA', 'dRNA']
    for i in range(3):
        data = [['Result'] + ngs]
        for j in range(len(r)):
            t = r[j]
            r1, r2 = compare(ngs, t[i], result)
            data.append([ls[j]] + r1)
            data.append([ls[j]] + r2)

        workbook = openpyxl.Workbook()
        sheet = workbook.active
        
        for col, column_data in enumerate(data, start=1):
            for row, value in enumerate(column_data, start=1):
                sheet.cell(row=row, column=col, value=value)
        workbook.save(outputdir + 'ngs_MCF7_' + list_types[i] + '.xlsx')