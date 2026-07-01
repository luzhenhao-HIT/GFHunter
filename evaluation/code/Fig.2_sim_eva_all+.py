# Calculate recall, precision, and F1 score of simulation data.
# All tools retains fusions with supporting reads number >= 3, 5, 8, 12 in 10x, 20x, 30x, 50x data respectively.
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
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            gene1 = parts[0]
            gene2 = parts[1]
            r = Result(gene1, gene2)
            key = create_key(gene1, gene2)
            result[key] = r
        line = file.readline()
    return result

def read_GFHunter(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/GFHunter.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            g1 = parts[0]
            g2 = parts[1]
            key = create_key(g1, g2)
            num = int(parts[7])
            if num >= n:
                if 'flag = 5' in line:
                    if key not in genefusion_exact.keys():
                        genefusion_exact[key] = GeneFuison(g1, g2, num)

                if 'flag = 5' in line or 'flag = 4' in line:
                    if key not in genefusion_approximate.keys():
                        genefusion_approximate[key] = GeneFuison(g1, g2, num)

        line = file.readline()
    return (genefusion_approximate, genefusion_exact)

def read_LongGF(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/LongGF.log')
    line = file.readline()
    longgf = {}
    while line:
        if 'SumGF' in line:
            parts = line.split('\t')[1].split(' ')
            if len(parts) >= 2:
                g1 = parts[0].split(':')[0]
                g2 = parts[0].split(':')[1]
                num = int(parts[1])
                key = create_key(g1, g2)
                if num >= n:
                    longgf[key] = GeneFuison(g1, g2, num)
        line = file.readline()
    return longgf

def read_JAFFAL(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/JAFFAL.csv')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        parts = line.strip().split(',')
        if len(parts) > 10:
            if 'Confidence' in line:
                g1 = parts[1].split(':')[0]
                g2 = parts[1].split(':')[1]
                num = int(parts[10])
                key = create_key(g1, g2)
                if num >= n:
                    jaffal[key] = GeneFuison(g1, g2, num)
        line = file.readline()
    return jaffal

def read_fusionseeker(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/fusionseeker.txt')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) > 3:
            g1 = parts[1]
            g2 = parts[2]
            num = int(parts[3])
            key = create_key(g1, g2)
            if num >= n:
                jaffal[key] = GeneFuison(g1, g2, num)
        line = file.readline()
    return jaffal

def read_Genion(depth, tp, n):
    file = open(path + '/result/All-Simulation/' + tp + str(depth) + 'x/Genion.tsv')
    line = file.readline()
    jaffal = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) > 4:
            g1 = parts[1].split('::')[0]
            g2 = parts[1].split('::')[1]
            num = int(parts[4])
            key = create_key(g1, g2)
            if num >= n:
                jaffal[key] = GeneFuison(g1, g2, num)
        line = file.readline()
    return jaffal

def read_ctat(depth, tp, n):
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
            num = int(parts[1]) # num_LR
            key = create_key(g1, g2)
            if num >= n:
                ctat[key] = GeneFuison(g1, g2, num)
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

def GFHunter_eva(depth, tp, result, n):
    genefusion_f4, genefusion_f5 = read_GFHunter(depth, tp, n)
    eva_f4 = calculate(result, genefusion_f4)
    eva_f5 = calculate(result, genefusion_f5)
    return ('GFHunter\t' + eva_f4, 'GFHunter_f5\t' + eva_f5)

def LongGF_eva(depth, tp, result, n):
    genefusion = read_LongGF(depth, tp, n)
    eva = calculate(result, genefusion)
    return 'LongGF\t' + eva

def JAFFAL_eva(depth, tp, result, n):
    genefusion = read_JAFFAL(depth, tp, n)
    eva = calculate(result, genefusion)
    return 'JAFFAL\t' + eva

def fusionseeker_eva(depth, tp, result, n):
    genefusion = read_fusionseeker(depth, tp, n)
    eva = calculate(result, genefusion)
    return 'fusionseeker\t' + eva

def Genion_eva(depth, tp, result, n):
    genefusion = read_Genion(depth, tp, n)
    eva = calculate(result, genefusion)
    return 'Genion\t' + eva

def ctat_eva(depth, tp, result, n):
    genefusion = read_ctat(depth, tp, n)
    eva = calculate(result, genefusion)
    return 'ctat-LR-fusion\t' + eva

def evaluation():
    tps = ['ONT', 'PB']
    depths = [10, 20, 30, 50]
    outputdir = path + '/result/Output/Fig.2 e&f +/'
    
    # 定义不同测序深度的 num 支持度阈值
    d = {10 : 3, 20 : 5, 30 : 8, 50 : 12}
    
    for tp in tps:
        for depth in depths:
            # 获取当前深度对应的支持reads过滤阈值 n，默认为 3
            n = d.get(depth, 3) 
            
            result = read_result(depth)
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)
                
            with open(outputdir + 'Sim_GF_'  + tp + str(depth) + 'x.tsv', 'w') as f:
                f.write('tools\tTP\tFP\tFN\trecall\tprecision\tF1-score\n')
                
                gf_f4, gf_f5 = GFHunter_eva(depth, tp, result, n)
                f.write(gf_f4)
                #f.write(gf_f5)
                
                lg = LongGF_eva(depth, tp, result, n)
                f.write(lg)
                
                ja = JAFFAL_eva(depth, tp, result, n)
                f.write(ja)
                
                fs = fusionseeker_eva(depth, tp, result, n)
                f.write(fs)
                
                ge = Genion_eva(depth, tp, result, n)
                f.write(ge)
                
                ct = ctat_eva(depth, tp, result, n)
                f.write(ct)

if __name__ == "__main__":
    evaluation()