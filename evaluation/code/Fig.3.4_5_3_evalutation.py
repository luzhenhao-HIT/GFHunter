# Calculate recall, precision, F1-score in real cell lines data: HCT-116 dRNA; SKBR-3 Iso-seq; MCF-7 cDNA, dRNA, Iso-seq when setting "Shared" fusions as ground truth.

import os
from collections import defaultdict
path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.3 e/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

result_dict = defaultdict(list)
t_dict = defaultdict(int)
all_dict = defaultdict(int)

def Genion(d, tp):
    dir = path + '/result/' + str(d) + '/'
    resultfile = open(dir + 'Genion_' + tp + '.tsv')
    global result_dict
    d = defaultdict(int)
    line = resultfile.readline()
    while line:
        gene1 = line.split('\t')[1].split('::')[0]
        gene2 = line.split('\t')[1].split('::')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['Genion'] += 1
        result_dict[key].append('Genion')

def JAFFAL_C(d, tp):
    dir = path + '/result/' + str(d) + '/'
    resultfile = open(dir + 'JAFFAL_' + tp + '.csv')
    global result_dict
    d = defaultdict(int)
    line = resultfile.readline()
    line = resultfile.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['JAFFAL_C'] += 1
        result_dict[key].append('JAFFAL_C')

def LongGF(d, tp):
    dir = path + '/result/' + str(d) + '/'
    resultfile = open(dir + 'LongGF_' + tp + '.log')
    global result_dict
    d = defaultdict(int)
    line = resultfile.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['LongGF'] += 1
        result_dict[key].append('LongGF')

def fusionseeker(d, tp):
    dir = path + '/result/' + str(d) + '/'
    resultfile = open(dir + 'fusionseeker_' + tp + '.txt')
    global result_dict
    d = defaultdict(int)
    line = resultfile.readline()
    line = resultfile.readline()
    while line:
        gene1 = line.split('\t')[1]
        gene2 = line.split('\t')[2]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['fusionseeker'] += 1
        result_dict[key].append('fusionseeker')

def GFHunter(d, tp):
    dir = path + '/result/' + str(d) + '/'
    resultfile = open(dir + 'GFHunter_' + tp + '.csv')
    global result_dict
    d = defaultdict(int)
    line = resultfile.readline()
    line = resultfile.readline()
    while line:
        if 'flag = 5' in line or 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['GFHunter_RF'] += 1
        result_dict[key].append('GFHunter_RF')

def compare(tresult):
    global result_dict, t_dict, all_dict
    TN = 0
    with open(tresult, 'w') as f:
        for key, value in result_dict.items():
            line = key
            if len(value) >= 3:
                for i in value:
                    t_dict[i] += 1
                    line += ',' + i
                f.write(line + '\n')
                TN += 1
        for key in ['GFHunter_RF', 'LongGF', 'JAFFAL_C', 'fusionseeker', 'Genion']:
            if key in t_dict.keys():
                TP = t_dict[key]
            else:
                TP = 0
            NP = all_dict[key]
            FP = NP - TP
            recall = float(TP/(TP+FP))
            precress = float(TP/TN)
            F1 = 2*recall*precress / (recall +precress)
            f.write(key + ',' + str(TP) + ',' + str(FP) + ',' + str(TN-TP) + ',' + str(recall) + ',' + str(precress) + ',' + str(F1)+ '\n')
    result_dict = defaultdict(list)
    t_dict = defaultdict(int)
    all_dict = defaultdict(int)

if __name__ == "__main__":
    dict = {'HCT-116':['dRNA'], 'SKBR-3':['PacBio'], 'MCF-7':['cDNA', 'dRNA', 'PacBio'] }
    for key, tps in dict.items():
        dir = str(key)
        for tp in tps:
            Genion(dir, tp)
            JAFFAL_C(dir, tp)
            LongGF(dir, tp)
            fusionseeker(dir, tp)
            GFHunter(dir, tp)
            compare(outputdir + str(key) + ' ' + tp + '.csv')