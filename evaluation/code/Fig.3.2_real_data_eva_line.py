# Calculate fusions number in real cell lines data: HCT-116 cDNA, dRNA; SKBR-3 Iso-seq; MCF-7 cDNA, dRNA, Iso-seq.
# GFHunter retains fusions with supporting reads number >= 4 in the large datasets HCT-116 cDNA.
# Note1: GFHunter's results are default supporting reads number >= 2, so we sort the reads number in ths program. 
# Note2: When using GFHunter to detect fusions, users can change "-l num(int)" to limit the supporting reads number.
# Note3: In MCF-7 cDNA and dRNA datasets, GFHunter can find one more fusions validated when "--max_exon_boundary 50" is setted, result of which marked by "_e50".
import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def read_GFHunter(d, tp):
    dir = path + '/result/' + str(d) + '/'
    if d == 'MCF-7' and (tp == 'cDNA' or tp == 'dRNA'):
        file = open(dir + 'GFHunter_' + tp + '_e50.csv')
    else:
        file = open(dir + 'GFHunter_' + tp + '.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    genefusion_all = {}
    while line:
        g1 = line.split('\t')[0]
        g2 = line.split('\t')[1]
        num = int(line.split('\t')[7])
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        if d == 'HCT-116' and tp == 'cDNA':
            if num >= 4:
                if 'flag = 5' in line:
                    genefusion_exact[key] = 1
                if 'flag = 4' in line or 'flag = 5' in line:
                    genefusion_approximate[key] = 1
                if 'flag = 3' in line or 'flag = 4' in line or 'flag = 5' in line:
                    genefusion_all[key] = 1
        else:
            if 'flag = 5' in line:
                genefusion_exact[key] = 1
            if 'flag = 4' in line or 'flag = 5' in line:
                genefusion_approximate[key] = 1
            if 'flag = 3' in line or 'flag = 4' in line or 'flag = 5' in line:
                genefusion_all[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    l2 = len(genefusion_approximate.keys()) - l1
    l3 = len(genefusion_all.keys()) -l1 - l2
    print('GFHunter: ' ,l1, l2, l3)

def read_LongGF(d, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'LongGF_' + tp + '.log')
    line = file.readline()
    genefusion_exact = {}
    while line:
        if 'SumGF' in line:
            g1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            g2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            l = [g1, g2]
            l.sort()
            g1 = l[0]
            g2 = l[1]
            key = g1 + ':' + g2
            genefusion_exact[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    print('LongGF: ', l1)

def read_JAFFAL(d, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'JAFFAL_' + tp + '.csv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    genefusion_approximate = {}
    genefusion_all = {}
    while line:
        g1 = line.split(',')[1].split(':')[0]
        g2 = line.split(',')[1].split(':')[1]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        if 'HighConfidence' in line:
            genefusion_exact[key] = 1
        if 'LowConfidence' in line or 'HighConfidence' in line:
            genefusion_approximate[key] = 1
        if 'PotentialTransSplicing' in line or 'LowConfidence' in line or 'HighConfidence' in line:
            genefusion_all[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    l2 = len(genefusion_approximate.keys()) - l1
    l3 = len(genefusion_all.keys()) - l2 - l1
    print('JAFFAL: ', l1, l2, l3)

def read_FusionSeeker(d, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'fusionseeker_' + tp + '.txt')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    while line:
        g1 = line.split('\t')[1]
        g2 = line.split('\t')[2]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        genefusion_exact[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    print('FusiomSeeker: ', l1)

def read_Genion(d, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'Genion_' + tp + '.tsv')
    line = file.readline()
    line = file.readline()
    genefusion_exact = {}
    while line:
        g1 = line.split('\t')[1].split('::')[0]
        g2 = line.split('\t')[1].split('::')[1]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        genefusion_exact[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    print('Genion: ', l1)

# 新增：读取并统计 ctat-LR-fusion 中基因对总数
def read_ctat(d, tp):
    dir = path + '/result/' + str(d) + '/'
    try:
        file = open(dir + 'ctat-LR-fusion_' + tp + '.tsv')
    except FileNotFoundError:
        try:
            file = open(dir + 'ctat-LR-fusion.' + tp + '.tsv')
        except FileNotFoundError:
            file = open(dir + 'ctat-LR-fusion.tsv')
            
    line = file.readline() # 跳过表头
    line = file.readline()
    genefusion_exact = {}
    while line:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            g1 = parts[2]
            g2 = parts[5]
            l = [g1, g2]
            l.sort()
            g1 = l[0]
            g2 = l[1]
            key = g1 + ':' + g2
            genefusion_exact[key] = 1
        line = file.readline()
    l1 = len(genefusion_exact.keys())
    print('ctat-LR-fusion: ', l1)

if __name__ == "__main__":
    dict = {'HCT-116':['cDNA', 'dRNA'], 'SKBR-3':['PacBio'], 'MCF-7':['cDNA', 'dRNA', 'PacBio'] }
    for key, tps in dict.items():
        dir = str(key)
        for tp in tps:
            print(dir, tp)
            read_GFHunter(dir, tp)
            read_LongGF(dir, tp)
            read_JAFFAL(dir, tp)
            read_FusionSeeker(dir, tp)
            read_Genion(dir, tp)
            # 调用 ctat 的解析方法打印总数
            read_ctat(dir, tp)