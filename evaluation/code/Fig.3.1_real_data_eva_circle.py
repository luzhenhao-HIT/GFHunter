# Calculate validated fusions in real cell lines data: HCT-116 cDNA, dRNA; SKBR-3 Iso-seq; MCF-7 cDNA, dRNA, Iso-seq.
# GFHunter retains fusions with supporting reads number >= 4 in the large datasets HCT-116 cDNA.
# Note1: GFHunter's results are default supporting reads number >= 2, so we sort the reads number in ths program. 
# Note2: When using GFHunter to detect fusions, users can change "-l num(int)" to limit the supporting reads number.
# Note3: In MCF-7 cDNA and dRNA datasets, GFHunter can find one more fusions validated when "--max_exon_boundary 50" is setted, result of which marked by "_e50".
import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.3 c1/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

class Result():
    def __init__(self, g1, g2, name):
        self.g1 = g1
        self.g2 = g2
        self.name = name
        pass

def read_result(d):
    dir = path + '/result/' + str(d) + '/'
    results = {}
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
        r = Result(g1, g2, name)
        results[key] = r
        line = file.readline()
    return results

def read_GFHunter(d, results, tp):
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
                    if key in results.keys():
                        genefusion_exact[key] = 1
                    else:
                        genefusion_exact[key] = 0

                if 'flag = 5' in line or 'flag = 4' in line:
                    if key in results.keys():
                        genefusion_approximate[key] = 1
                    else:
                        genefusion_approximate[key] = 0
                if key in results.keys():
                    genefusion_all[key] = 1
                else:
                    genefusion_all[key] = 0
        else:
            if 'flag = 5' in line:
                if key in results.keys():
                    genefusion_exact[key] = 1
                else:
                    genefusion_exact[key] = 0

            if 'flag = 5' in line or 'flag = 4' in line:
                if key in results.keys():
                    genefusion_approximate[key] = 1
                else:
                    genefusion_approximate[key] = 0
            if key in results.keys():
                genefusion_all[key] = 1
            else:
                genefusion_all[key] = 0

        line = file.readline()
    with open(outputdir + d + ' ' + tp + ' GFHunter RF.csv', 'w') as f:
        f.write('Fusion\tGFHunter\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_exact.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))
    with open(outputdir + d + ' ' + tp + ' GFHunter RF+SF.csv', 'w') as f:
        f.write('Fusion\tGFHunter\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_approximate.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))
    with open(outputdir + d + ' ' + tp + ' GFHunter RF+SF+PF.csv', 'w') as f:
        f.write('Fusion\tGFHunter\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_all.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))

def read_LongGF(d, results, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'LongGF_' + tp + '.log')
    line = file.readline()
    genefusion = {}
    while line:
        if 'SumGF' in line:
            g1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            g2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            l = [g1, g2]
            l.sort()
            g1 = l[0]
            g2 = l[1]
            key = g1 + ':' + g2
            if key in results.keys():
                genefusion[key] = 1
            else:
                genefusion[key] = 0

        line = file.readline()
    with open(outputdir + d + ' ' + tp + ' LongGF.csv', 'w') as f:
        f.write('Fusion\tLongGF\n')
        n = 0
        for key, item in results.items():
            if key in genefusion.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))

def read_JAFFAL(d, results, tp):
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
            if key in results.keys():
                genefusion_exact[key] = 1
            else:
                genefusion_exact[key] = 0
        if 'Confidence' in line:
            if key in results.keys():
                genefusion_approximate[key] = 1
            else:
                genefusion_approximate[key] = 0

        if key in results.keys():
            genefusion_all[key] = 1
        else:
            genefusion_all[key] = 0

        line = file.readline()
    with open(outputdir + d + ' ' + tp + ' JAFFAL HC.csv', 'w') as f:
        f.write('Fusion\tJAFFAL\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_exact.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))
    with open(outputdir + d + ' ' + tp + ' JAFFAL HC+LC.csv', 'w') as f:
        f.write('Fusion\tJAFFAL\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_approximate.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))
    with open(outputdir + d + ' ' + tp + ' JAFFAL HC+LC+PT.csv', 'w') as f:
        f.write('Fusion\tJAFFAL\n')
        n = 0
        for key, item in results.items():
            if key in genefusion_all.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))

def read_fusionseeker(d, results, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'fusionseeker_' + tp + '.txt')
    line = file.readline()
    line = file.readline()
    genefusion = {}
    while line:
        g1 = line.split('\t')[1]
        g2 = line.split('\t')[2]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        if key in results.keys():
            genefusion[key] = 1
        else:
            genefusion[key] = 0

        line = file.readline()
    with open(outputdir + d + ' ' + tp + ' fusionseeker.csv', 'w') as f:
        f.write('Fusion\tfusionseeker\n')
        n = 0
        for key, item in results.items():
            if key in genefusion.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))

def read_genion(d, results, tp):
    dir = path + '/result/' + str(d) + '/'
    file = open(dir + 'Genion_' + tp + '.tsv')
    line = file.readline()
    line = file.readline()
    genefusion = {}
    while line:
        g1 = line.split('\t')[1].split('::')[0]
        g2 = line.split('\t')[1].split('::')[1]
        l = [g1, g2]
        l.sort()
        g1 = l[0]
        g2 = l[1]
        key = g1 + ':' + g2
        if key in results.keys():
            genefusion[key] = 1
        else:
            genefusion[key] = 0

        line = file.readline()
    with open(outputdir + d + ' ' + tp + ' Genion.csv', 'w') as f:
        f.write('Fusion\tGenion\n')
        n = 0
        for key, item in results.items():
            if key in genefusion.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))

# 新增：读取 ctat-LR-fusion，结构与前述一致
def read_ctat(d, results, tp):
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
    genefusion = {}
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
            
            if key in results.keys():
                genefusion[key] = 1
            else:
                genefusion[key] = 0
        line = file.readline()
        
    with open(outputdir + d + ' ' + tp + ' ctat.csv', 'w') as f:
        f.write('Fusion\tctat-LR-fusion\n')
        n = 0
        for key, item in results.items():
            if key in genefusion.keys():
                f.write(item.name + '\t1\n')
                n += 1
            else:
                f.write(item.name + '\n')
        f.write('total\t' + str(n))


if __name__ == "__main__":
    dict = {'HCT-116':['cDNA', 'dRNA'], 'SKBR-3':['PacBio'], 'MCF-7':['cDNA', 'dRNA', 'PacBio'] }
    for key, ls in dict.items():
        dir = str(key)
        results = read_result(dir)
        for tp in ls:
            read_GFHunter(dir, results, tp)
            read_LongGF(dir, results, tp)
            read_JAFFAL(dir, results, tp)
            read_fusionseeker(dir, results, tp)
            read_genion(dir, results, tp)
            # 增加对 ctat 评估调用的代码
            read_ctat(dir, results, tp)