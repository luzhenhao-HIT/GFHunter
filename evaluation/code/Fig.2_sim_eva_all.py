# Calculate recall, precision, and F1 score of simulation data.
# GFHunter retains fusions with supporting reads number >= 3, 5, 8, 12 in 10x, 20x, 30x, 50x data respectively.
# Note1: GFHunter's results are default supporting reads number >= 2, so we sort the reads number in ths program. 
# Note2: When using GFHunter to detect fusions, users can change "-l num(int)" to limit the supporting reads number.

import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class Result():
    def __init__(self, g1, g2):
        self.g1 = g1
        self.g2 = g2
        self.key = False
        pass

class GeneFuison():
    def __init__(self, g1, g2, num):
        self.g1 = g1
        self.g2 = g2
        self.num = num
        self.key = False
        self.bpkey = False
        pass

def create_key(a, b):
    ls = [a, b]
    ls.sort()
    key = ls[0] + ':' + ls[1]
    return key

def read_result(depth):
    file = open(path + '/result/All-Simulation/simulation_result.tsv')
    line = file.readline()
    line = file.readline()
    result = {}
    while line:
        gene1 = line.split('\t')[0]
        gene2 = line.split('\t')[1]
        line = file.readline()
        r = Result(gene1, gene2)
        key = create_key(gene1, gene2)
        result[key] = r
    return result

def read_GFHunter(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/GFHunter.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    while line:
        g1 = line.split('\t')[0]
        g2 = line.split('\t')[1]
        key = create_key(g1, g2)
        num = int(line.split('\t')[7])
        if num >= n:
            if 'flag = 5' in line:
                if key not in genefusion_exact.keys():
                    genefusion_exact[key] = GeneFuison(g1, g2, num)

            if 'flag = 5' in line or 'flag = 4' in line:
                if key not in genefusion_approximate.keys():
                    genefusion_approximate[key] = GeneFuison(g1, g2, num)

        line = file.readline()
    return (genefusion_approximate, genefusion_exact)

def read_LongGF(depth, tp):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/LongGF.log')
    line = file.readline()
    longgf = {}
    while line:
        if 'SumGF' in line:
            g1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            g2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            key = create_key(g1, g2)
            longgf[key] = GeneFuison(g1, g2, 0)
        line = file.readline()
    return longgf

def read_JAFFAL(depth, tp):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/JAFFAL.csv')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        if 'Confidence' in line:
            g1 = line.split(',')[1].split(':')[0]
            g2 = line.split(',')[1].split(':')[1]
            key = create_key(g1, g2)
            jaffal[key] = GeneFuison(g1, g2, 0)
        line = file.readline()
    return jaffal

def read_fusionseeker(depth, tp):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/fusionseeker.txt')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1]
        g2 = line.split('\t')[2]
        key = create_key(g1, g2)
        jaffal[key] = GeneFuison(g1, g2, 0)
        line = file.readline()
    return jaffal

def read_Genion(depth, tp):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/Genion.tsv')
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1].split('::')[0]
        g2 = line.split('\t')[1].split('::')[1]
        key = create_key(g1, g2)
        jaffal[key] = GeneFuison(g1, g2, 0)
        line = file.readline()
    return jaffal

# 新增：读取 ctat-LR-fusion 格式文件
def read_ctat(depth, tp):
    # 使用 try-except 兼容不同可能的文件命名习惯
    try:
        file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/ctat-LR-fusion.tsv')
    except FileNotFoundError:
        try:
            file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/ctat-LR-fusion_' + tp + str(depth) + 'x.tsv')
        except FileNotFoundError:
            file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/ctat-LR-fusion.' + tp + str(depth) + 'x.tsv')
            
    line = file.readline() # 跳过 header
    line = file.readline()
    ctat = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            g1 = parts[2] # LeftGene
            g2 = parts[5] # RightGene
            key = create_key(g1, g2)
            ctat[key] = GeneFuison(g1, g2, 0)
        line = file.readline()
    return ctat

def calculate(result, genefusion):
    TP = 0
    FP = 0
    FN = 0
    for key, res in result.items():
        if key in genefusion.keys():
            res.key = True
            genefusion[key].key = True
            
    for key, item in genefusion.items():
        if item.key == True:
            TP += 1
        else:
            FP += 1
            
    for key, item in result.items():
        if item.key == False:
            FN += 1
            
    # 防止分母为 0 报错 (ZeroDivisionError)
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    f1 = 2 * recall * precision / (recall + precision) if (recall + precision) > 0 else 0.0
    
    for key, item in genefusion.items():
        item.key = False
        item.bpkey = False
    for key, item in result.items():
        item.key = False
        
    line = str(TP) + '\t' + str(FP) + '\t' + str(FN) + '\t' + str(recall) + '\t' + str(precision) + '\t' + str(f1) + '\n'
    return line

def GFHunter_eva(depth, tp, result):
    d = {10 : 3, 20 : 5, 30 : 8, 50 : 12}
    genefusion_f4, genefusion_f5 = read_GFHunter(depth, tp, d[depth])
    eva_f4 = calculate(result, genefusion_f4)
    eva_f5 = calculate(result, genefusion_f5)
    return ('GFHunter\t' + eva_f4, 'GFHunter_f5\t' + eva_f5)

def LongGF_eva(depth, tp, result):
    genefusion = read_LongGF(depth, tp)
    eva = calculate(result, genefusion)
    return 'LongGF\t' + eva

def JAFFAL_eva(depth, tp, result):
    genefusion = read_JAFFAL(depth, tp)
    eva = calculate(result, genefusion)
    return 'JAFFAL\t' + eva

def fusionseeker_eva(depth, tp, result):
    genefusion = read_fusionseeker(depth, tp)
    eva = calculate(result, genefusion)
    return 'fusionseeker\t' + eva

def Genion_eva(depth, tp, result):
    genefusion = read_Genion(depth, tp)
    eva = calculate(result, genefusion)
    return 'Genion\t' + eva

# 新增：封装 ctat-LR-fusion 评估逻辑
def ctat_eva(depth, tp, result):
    genefusion = read_ctat(depth, tp)
    eva = calculate(result, genefusion)
    return 'ctat-LR-fusion\t' + eva

def evaluation():
    tps = ['ONT', 'PB']
    depths = [10, 20, 30, 50]
    outputdir = path + '/result/Output/Fig.2 e&f/'
    for tp in tps:
        for depth in depths:
            result = read_result(depth)
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)
                
            # 这里将 .csv 结尾修改为了 .tsv 结尾
            with open(outputdir + 'Sim_GF_'  + tp + str(depth) + 'x.tsv', 'w') as f:
                # 修复了 NF 拼写错误为 FN 
                f.write('tools\tTP\tFP\tFN\trecall\tprecision\tF1-score\n')
                
                gf_f4, gf_f5 = GFHunter_eva(depth, tp, result)
                f.write(gf_f4)
                #f.write(gf_f5)
                
                lg = LongGF_eva(depth, tp, result)
                f.write(lg)
                
                ja = JAFFAL_eva(depth, tp, result)
                f.write(ja)
                
                fs = fusionseeker_eva(depth, tp, result)
                f.write(fs)
                
                ge = Genion_eva(depth, tp, result)
                f.write(ge)
                
                # 写出 ctat 评估结果
                ct = ctat_eva(depth, tp, result)
                f.write(ct)

if __name__ == "__main__":
    evaluation()