# Calculate the consistency in real cell lines data: HCT-116.

import os
import openpyxl

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.4 b/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)
 
def GFHunter_RF(result_dir):
    l2 = []
    l3 = []
    cDNA_file = open(result_dir + 'GFHunter_cDNA.csv')
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def GFHunter_SF(result_dir):
    l2 = []
    l3 = []
    cDNA_file = open(result_dir + 'GFHunter_cDNA.csv')
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def GFHunter(result_dir):
    l2 = []
    l3 = []
    cDNA_file = open(result_dir + 'GFHunter_cDNA.csv')
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def Genion(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def JAFFAL_HC(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def JAFFAL_C(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data


def JAFFAL(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def LongGF(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def fusionseeker(result_dir):
    l2 = []
    l3 = []
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
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l2.sort()
    l3.sort()
    data = (l2, l3)
    print(len(l2), len(l3))
    return data

def JAFFA():
    l1 = []
    PacBio_file = open(path + '/result/NGS/JAFFA_Hct116.csv')
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
    PacBio_file = open(path + '/result/NGS/starfusion_Hct116.tsv')
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
    dir = path + '/result/HCT-116/'
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
    r = (gf_rf, gf_sf, gf, lg, ja_hc, ja_c, ja, fs, ge)
    jan = JAFFA()
    sf = starfusion()
    ngs = []
    for key in jan:
        if key in sf:
            ngs.append(key)
    ls = ['GFHunter RF', 'GFHunter RF + SF', 'GFHunter', 'LongGF', 'JAFFAL_HC', 'JAFFAL_C', 'JAFFAL', 'fusionseeker', 'Genion']
    list = ['cDNA', 'dRNA']
    for i in range(2):
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
        workbook.save(outputdir + 'ngs_HCT116_' + list[i] + '.xlsx')