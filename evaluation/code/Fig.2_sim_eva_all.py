# Calculate recall, precision, and F1 score of simulation data.
# GFHunter retains fusions with supporting reads number >= 3, 5, 8, 12 in 10x, 20x, 30x, 50x data respectively.
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
    file = open(path + '/result/All-Simualtion/result' + str(depth) + 'x.txt')
    line = file.readline()
    result = {}
    while line:
        if '>' in line:
            line = file.readline()
            gene1 = line.split('\t')[1]
            breakpoint1 = int(line.split('\t')[9])
            line = file.readline()
            gene2 = line.split('\t')[1]
            breakpoint2 = int(line.split('\t')[9])
            line = file.readline()
            r = Result(gene1, gene2, breakpoint1, breakpoint2)
            key = create_key(gene1, gene2)
            result[key] = r
    return result

def read_GFHunter(depth, tp, n):
    file = open(path + '/result/All-Simualtion/' + tp + str(depth) + 'x/GFHunter.csv')
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

def read_LongGF(depth, tp):
    file = open(path + '/result/All-Simualtion/' + tp + str(depth) + 'x/LongGF.log')
    line = file.readline()
    longgf = {}
    while line:
        if 'SumGF' in line:
            g1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            g2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            bp1 = int(line.split('\t')[1].split(' ')[2].split(':')[1])
            bp2 = int(line.split('\t')[1].split(' ')[3].split('\n')[0].split(':')[1])
            key = create_key(g1, g2)
            longgf[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return longgf

def read_JAFFAL(depth, tp):
    file = open(path + '/result/All-Simualtion/' + tp + str(depth) + 'x/JAFFAL.csv')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        if 'Confidence' in line:
            g1 = line.split(',')[1].split(':')[0]
            g2 = line.split(',')[1].split(':')[1]
            bp1 = int(line.split(',')[3])
            bp2 = int(line.split(',')[6])
            key = create_key(g1, g2)
            #print(key)
            jaffal[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return jaffal

def read_fusionseeker(depth, tp):
    file = open(path + '/result/All-Simualtion/' + tp + str(depth) + 'x/fusionseeker.txt')
    line = file.readline()
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1]
        g2 = line.split('\t')[2]
        bp1 = int(line.split('\t')[5])
        bp2 = int(line.split('\t')[7])
        key = create_key(g1, g2)
        #print(key)
        jaffal[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return jaffal

def read_Genion(depth, tp):
    file = open(path + '/result/All-Simualtion/' + tp + str(depth) + 'x/Genion.tsv')
    line = file.readline()
    jaffal = {}
    while line:
        g1 = line.split('\t')[1].split('::')[0]
        g2 = line.split('\t')[1].split('::')[1]
        bp1 = 0
        bp2 = 0
        key = create_key(g1, g2)
        #print(key)
        jaffal[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return jaffal

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
    #print(TP, FP, FN)
    recall = TP/(TP + FN)
    precision = TP/(TP + FP)
    f1 = 2 * recall * precision/(recall + precision)
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

def evaluation():
    tps = ['ONT', 'PB']
    depths = [10, 20, 30, 50]
    outputdir = path + '/result/Output/Fig.2 e&f/'
    for tp in tps:
        for depth in depths:
            result = read_result(depth)
            if not os.path.exists(outputdir):
                os.makedirs(outputdir)
            with open(outputdir + 'Sim_GF_'  + tp + str(depth) + 'x.csv', 'w') as f:
                f.write('tools\tTP\tFP\tNF\trecall\tprecision\tF1-score\n')
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
            


if __name__ == "__main__":
    evaluation()