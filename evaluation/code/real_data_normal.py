# Calculate false positive number of non-tumor real datasets.
# All tools retains fusions with supporting reads number >= 3, 8, 3 in Iso-seq, cDNA, dRNA data.
# Note1: All tools results are default supporting reads number >= 2, so we sort the reads number in ths program. 

import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Table 1/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

K = 50
d = {'PB' : 5, 'cDNA' : 8, 'dRNA' : 3}

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

def read_result():
    result = {}
    return result

def read_GFHunter(tp, n):
    file = open(path  + '/result/Non-tumor-data/GFHunter_' + tp + '.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    genefusion_all = {}
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
            
            if key in genefusion_all.keys():
                genefusion_all[key].add_bp(bp1, bp2)
            else:
                genefusion_all[key] = GeneFuison(g1, g2, bp1, bp2, num)

        line = file.readline()
    return (genefusion_approximate, genefusion_exact, genefusion_all)

def read_LongGF(tp, n):
    file = open(path  + '/result/Non-tumor-data/LongGF_' + tp + '.log')
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
                longgf[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return longgf

def read_JAFFAL(tp, n):
    file = open(path  + '/result/Non-tumor-data/JAFFAL_' + tp + '.csv')
    line = file.readline()
    line = file.readline()
    jaffal_HC = {}
    jaffal_LC = {}
    jaffal_PT = {}
    while line:
        g1 = line.split(',')[1].split(':')[0]
        g2 = line.split(',')[1].split(':')[1]
        num = int(line.split(',')[10])
        bp1 = int(line.split(',')[3])
        bp2 = int(line.split(',')[6])
        key = create_key(g1, g2)
        if num >= n:
            if 'HighConfidence' in line:
                jaffal_HC[key] = GeneFuison(g1, g2, bp1, bp2, 0)
            if 'HighConfidence' in line or 'LowConfidence' in line:
                jaffal_LC[key] = GeneFuison(g1, g2, bp1, bp2, 0)
            jaffal_PT[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return (jaffal_PT, jaffal_LC, jaffal_HC)

def read_fusionseeker(tp, n):
    file = open(path  + '/result/Non-tumor-data/fusionseeker_' + tp + '.txt')
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
            jaffal[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return jaffal

def read_Genion(tp, n):
    file = open(path  + '/result/Non-tumor-data/Genion_' + tp + '.tsv')
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
            jaffal[key] = GeneFuison(g1, g2, bp1, bp2, 0)
        line = file.readline()
    return jaffal

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
    line = str(FP) + '\n'
    return (line, keys)

def GFHunter_eva(tp, result):
    genefusion_f4, genefusion_f5, genefusion_f3 = read_GFHunter(tp, d[tp])
    eva_f4, keys_4 = calculate(result, genefusion_f4)
    eva_f5, keys_5 = calculate(result, genefusion_f5)
    eva_f3, keys_3 = calculate(result, genefusion_f3)
    return ('GFHunter_SF\t' + eva_f4, keys_4, 'GFHunter_RF\t' + eva_f5, keys_5, 'GFHunter_PF\t' + eva_f3, keys_3)

def LongGF_eva(tp, result):
    genefusion = read_LongGF(tp, d[tp])
    eva, keys = calculate(result, genefusion)
    return ('LongGF\t' + eva, keys)

def JAFFAL_eva( tp, result):
    jaffal_PT, jaffal_LC, jaffal_HC = read_JAFFAL(tp, d[tp])
    eva_PT, keys_PT = calculate(result, jaffal_PT)
    eva_LC, keys_LC = calculate(result, jaffal_LC)
    eva_HC, keys_HC = calculate(result, jaffal_HC)
    return ('JAFFAL_HC\t' + eva_HC, keys_HC, 'JAFFAL_LC\t' + eva_LC, keys_LC, 'JAFFAL_PT\t' + eva_PT, keys_PT)

def fusionseeker_eva(tp, result):
    genefusion = read_fusionseeker(tp, d[tp])
    eva, keys = calculate(result, genefusion)
    return ('fusionseeker\t' + eva, keys)

def Genion_eva(tp, result):
    genefusion = read_Genion(tp, d[tp])
    eva, keys = calculate(result, genefusion)
    return ('Genion\t' + eva, keys)


def evaluation():
    tps = ['PB', 'cDNA', 'dRNA']
    GF = []
    LG = []
    JA = []
    FS = []
    GE = []
    for tp in tps:
        result = read_result()
        with open(outputdir + 'non-tumor ' + tp + '.csv', 'w') as f:
            f.write('tools\tFP\n')
            gf_f4, k14, gf_f5, k1 ,gf_f3, k13 = GFHunter_eva(tp, result)
            f.write(gf_f5)
            f.write(gf_f4)
            #f.write(gf_f3)
            lg, k2 = LongGF_eva(tp, result)
            f.write(lg)
            ja_hc, k3, ja_lc, k32, ja_pt, k33 = JAFFAL_eva(tp, result)
            f.write(ja_hc)
            f.write(ja_lc)
            #f.write(ja_pt)
            fs, k4 = fusionseeker_eva(tp, result)
            f.write(fs)
            ge, k5 = Genion_eva(tp, result)
            f.write(ge)
            k0 = k1 + k2 + k3 + k4 + k5
            k = list(set(k0))
            keys = {}
            for key in k:
                n = []
                if key in k1:
                    n.append(1)
                if key in k2:
                    n.append(2)
                if key in k3:
                    n.append(3)
                if key in k4:
                    n.append(4)
                if key in k5:
                    n.append(5)
                keys[key] = n
            n = [0, 0, 0, 0, 0]
            for key, nums in keys.items():
                n[len(nums) - 1] += 1
            for i in n:
                f.write(str(i) + '\n')

if __name__ == "__main__":
    evaluation()