# Calculate false positive number of nagetive simulation data.
# GFHunter retains fusions with supporting reads number >= 8 in 30x data.
# Note1: GFHunter's results are default supporting reads number >= 2, so we sort the reads number in ths program. 
# Note2: When using GFHunter to detect fusions, users can change "-l num(int)" to limit the supporting reads number.

import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class Result():
    def __init__(self, g1, g2, bp1, bp2):
        self.g1 = g1
        self.g2 = g2
        self.bp1 = bp1
        self.bp2 = bp2
        self.key = False
        pass

class GeneFuison():
    def __init__(self, g1, g2, bp1, bp2, num):
        self.g1 = g1
        self.g2 = g2
        self.bp1s = [bp1,]
        self.bp2s = [bp2,]
        self.key = False
        self.bpkey = False
        self.num = num
        pass
    def add_bp(self, bp1, bp2):
        self.bp1s.append(bp1)
        self.bp2s.append(bp2)
        pass

def create_key(a, b):
    ls = [a, b]
    ls.sort()
    key = ls[0] + ':' + ls[1]
    return key

def read_result(depth):
    result = {}
    return result

def read_GFHunter(depth, tp, n):
    file = open(path + '/result/All-Simulation/Nagetive/GFHunter_' + tp + str(depth) + 'x.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    while line:
        g1 = line.split('\t')[0]
        g2 = line.split('\t')[1]
        bp1 = int(line.split('\t')[4].split('; ')[0].split(':')[1])
        bp2 = int(line.split('\t')[4].split('; ')[1].split(':')[1])
        key = create_key(g1, g2)
        num = int(line.split('\t')[7])
        
        if num >= n:
            if 'flag = 5' in line:
                if key in genefusion_exact.keys():
                    genefusion_exact[key].add_bp(bp1, bp2)
                else:
                    genefusion_exact[key] = GeneFuison(g1, g2, bp1, bp2, num)

            if 'flag = 5' in line or 'flag = 4' in line:
                if key in genefusion_approximate.keys():
                    genefusion_approximate[key].add_bp(bp1, bp2)
                else:
                    genefusion_approximate[key] = GeneFuison(g1, g2, bp1, bp2, num)

        line = file.readline()
    return (genefusion_approximate, genefusion_exact)

def read_LongGF(depth, tp, n):
    file = open(path + '/result/All-Simulation/Nagetive/LongGF_' + tp + str(depth) + 'x.log')
    line = file.readline()
    longgf = {}
    while line:
        if 'SumGF' in line:
            g1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            g2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            num = int(line.split('\t')[1].split(' ')[1])
            bp1 = int(line.split('\t')[1].split(' ')[2].split(':')[1])
            bp2 = int(line.split('\t')[1].split(' ')[3].split('\n')[0].split(':')[1])
            key = create_key(g1, g2)
            if num >= n:
                longgf[key] = GeneFuison(g1, g2, bp1, bp2, num)
        line = file.readline()
    return longgf

def read_JAFFAL(depth, tp, n):
    file = open(path + '/result/All-Simulation/Nagetive/JAFFAL_' + tp + str(depth) + 'x.csv') 
    line = file.readline()
    line = file.readline()
    jaffal = {}
    jaffal_HC = {}
    while line:
        parts = line.split(',')
        if len(parts) > 10:
            g1 = parts[1].split(':')[0]
            g2 = parts[1].split(':')[1]
            bp1 = int(parts[3])
            bp2 = int(parts[6])
            num = int(parts[10])
            key = create_key(g1, g2)
            
            if num >= n:
                if 'Confidence' in line:
                    jaffal[key] = GeneFuison(g1, g2, bp1, bp2, num)
                if 'HighConfidence' in line:
                    jaffal_HC[key] = GeneFuison(g1, g2, bp1, bp2, num)
        line = file.readline()
    return (jaffal, jaffal_HC)

def read_fusionseeker(depth, tp, n):
    file = open(path + '/result/All-Simulation/Nagetive/fusionseeker_' + tp + str(depth) + 'x.txt')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1]
        g2 = line.split('\t')[2]
        num = int(line.split('\t')[3])
        bp1 = int(line.split('\t')[5])
        bp2 = int(line.split('\t')[7])
        key = create_key(g1, g2)
        if num >= n:
            jaffal[key] = GeneFuison(g1, g2, bp1, bp2, num)
        line = file.readline()
    return jaffal

def read_Genion(depth, tp, n):
    file = open(path + '/result/All-Simulation/Nagetive/Genion_' + tp + str(depth) + 'x.tsv')
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1].split('::')[0]
        g2 = line.split('\t')[1].split('::')[1]
        num = int(line.split('\t')[4])
        bp1 = 0
        bp2 = 0
        key = create_key(g1, g2)
        if num >= n:
            jaffal[key] = GeneFuison(g1, g2, bp1, bp2, num)
        line = file.readline()
    return jaffal

def read_ctat(depth, tp, n):
    try:
        file = open(path + '/result/All-Simulation/Nagetive/ctat-LR-fusion.' + tp + str(depth) + 'x.tsv')
    except FileNotFoundError:
        file = open(path + '/result/All-Simulation/Nagetive/ctat-LR-fusion_' + tp + str(depth) + 'x.tsv')
        
    line = file.readline()
    line = file.readline()
    ctat = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            g1 = parts[2] 
            g2 = parts[5] 
            bp1 = int(parts[4].split(':')[1])
            bp2 = int(parts[7].split(':')[1]) 
            num = int(parts[1]) 
            key = create_key(g1, g2)
            
            if num >= n:
                ctat[key] = GeneFuison(g1, g2, bp1, bp2, num)
        line = file.readline()
    return ctat

def calculate(result, genefusion):
    TP = 0
    FP = 0
    FN = 0
    keys = []
    for key, res in result.items():
        if key in genefusion.keys():
            res.key = True
            genefusion[key].key = True
    for key, item in genefusion.items():
        if item.key == True:
            TP += 1
        else:
            FP += 1
            keys.append(key)
    for key, item in result.items():
        if item.key == False:
            FN += 1
    line = str(FP)
    return (line, keys)

def GFHunter_eva(depth, tp, result, n):
    genefusion_f4, genefusion_f5 = read_GFHunter(depth, tp, n)
    eva_f4, keys_4 = calculate(result, genefusion_f4)
    eva_f5, keys_5 = calculate(result, genefusion_f5)
    return ('GFHunter SF\t' + str(int(eva_f4) - int(eva_f5)), keys_4, 'GFHunter RF\t' + eva_f5, keys_5)

def LongGF_eva(depth, tp, result, n):
    genefusion = read_LongGF(depth, tp, n)
    eva, keys = calculate(result, genefusion)
    return ('LongGF\t' + eva, keys)

def JAFFAL_eva(depth, tp, result, n):
    genefusion_c, genefusion_hc = read_JAFFAL(depth, tp, n)
    eva, keys = calculate(result, genefusion_c)
    eva_h, keys_h = calculate(result, genefusion_hc)
    return ('JAFFAL HC+LC\t' + str(int(eva) - int(eva_h)), keys, 'JAFFAL HC\t' + eva_h, keys_h)

def fusionseeker_eva(depth, tp, result, n):
    genefusion = read_fusionseeker(depth, tp, n)
    eva, keys = calculate(result, genefusion)
    return ('fusionseeker\t' + eva, keys)

def Genion_eva(depth, tp, result, n):
    genefusion = read_Genion(depth, tp, n)
    eva, keys = calculate(result, genefusion)
    return ('Genion\t' + eva, keys)

def ctat_eva(depth, tp, result, n):
    genefusion = read_ctat(depth, tp, n)
    eva, keys = calculate(result, genefusion)
    return ('ctat-LR-fusion\t' + eva, keys)

def k(keys):
    if not keys: return ""
    k = keys[0]
    for key in keys[1:]:
        k += '; ' + key
    return k

def evaluation():
    tps = ['ONT', 'PB']
    depths = [30]
    outputdir = path + '/result/Output/Fig.2 g/'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    # 定义不同测序深度的 num 支持度阈值
    d = {10 : 3, 20 : 5, 30 : 8, 50 : 12}
    
    for tp in tps:
        for depth in depths:
            # 动态获取 n (支持 reads 数阈值)，如果depth不在d中默认取8
            n = d.get(depth, 8) 
            result = read_result(depth)
            
            with open(outputdir + 'Sim_Na_'  + tp + str(depth) + 'x.tsv', 'w') as f:
                f.write('tools\tFP\tFusions\n')
                gf_f4, k14, gf_f5, k15 = GFHunter_eva(depth, tp, result, n)
                f.write(gf_f5 + '\t' + k(k15) + '\n')
                f.write(gf_f4 + '\t' + k(list(set(k14) - set(k15))) + '\n')
                
                lg, k2 = LongGF_eva(depth, tp, result, n)
                f.write(lg + '\t' + k(k2) + '\n')
                
                ja, k3, ja2, k32 = JAFFAL_eva(depth, tp, result, n)
                f.write(ja2 + '\t' + k(k32) + '\n')
                f.write(ja + '\t' + k(list(set(k3) - set(k32))) + '\n')
                
                fs, k4 = fusionseeker_eva(depth, tp, result, n)
                f.write(fs + '\t' + k(k4) + '\n')
                
                ge, k5 = Genion_eva(depth, tp, result, n)
                f.write(ge + '\t' + k(k5) + '\n')
                
                ct, k6 = ctat_eva(depth, tp, result, n)
                f.write(ct + '\t' + k(k6) + '\n')

if __name__ == "__main__":
    evaluation()