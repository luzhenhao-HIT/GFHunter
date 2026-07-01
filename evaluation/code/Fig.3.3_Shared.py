# Calculate "Singleton", "Doubleton", "Shared" fusions number in real cell lines data: HCT-116 dRNA; SKBR-3 Iso-seq; MCF-7 cDNA, dRNA, Iso-seq.

import os
from collections import defaultdict
path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.3 d/'
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
        if 'HighConfidence' in line or 'LowConfidence' in line:
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
        if 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            d[key] += 1
        line = resultfile.readline()
    for key in d.keys():
        all_dict['GFHunter'] += 1
        result_dict[key].append('GFHunter')

def ctat(d, tp):
    dir = path + '/result/' + str(d) + '/'
    try:
        resultfile = open(dir + 'ctat-LR-fusion_' + tp + '.tsv')
    except FileNotFoundError:
        try:
            resultfile = open(dir + 'ctat-LR-fusion.' + tp + '.tsv')
        except FileNotFoundError:
            resultfile = open(dir + 'ctat-LR-fusion.tsv')
            
    global result_dict
    d_map = defaultdict(int)
    line = resultfile.readline() # header
    line = resultfile.readline()
    while line:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            gene1 = parts[2]
            gene2 = parts[5]
            num = int(parts[1]) # num_LR
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            if num >= 2:
                d_map[key] += 1
        line = resultfile.readline()
    for key in d_map.keys():
        all_dict['ctat'] += 1
        result_dict[key].append('ctat')

def compare(tresult):
    global result_dict, t_dict, all_dict
    shared = defaultdict(int)
    doubleton = defaultdict(int)
    singleton = defaultdict(int)
    with open(tresult, 'w') as f:
        for key, value in result_dict.items():
            # 修改为 >= 4 (6选4)
            if len(value) >= 4: 
                for i in value:
                    shared[i] += 1
            if len(value) == 2 or len(value) == 3:
                for i in value:
                    doubleton[i] += 1
            if len(value) == 1:
                for i in value:
                    singleton[i] += 1
        f.write('tools\tshared\trare\tsingleton\n')
        for key in ['GFHunter', 'LongGF', 'JAFFAL_C', 'fusionseeker', 'Genion', 'ctat']:
            line = key + '\t'
            if key in shared.keys():
                line += str(shared[key]) + '\t'
            else:
                line += '0\t'
            if key in doubleton.keys():
                line += str(doubleton[key]) + '\t'
            else:
                line += '0\t'
            if key in singleton.keys():
                line += str(singleton[key]) + '\n'
            else:
                line += '0\n'
            f.write(line)
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
            ctat(dir, tp) 
            compare(outputdir + str(key) + ' ' + tp + '.csv')