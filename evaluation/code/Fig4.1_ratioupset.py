# Calculate the counts and harmonic means of consistently reported fusions compared to total reported fusions in real cell lines data: MCF-7 cDNA, dRNA, Iso-seq.

import os

path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
outputdir = path + '/result/Output/Fig.4 a/'
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

def read_result():
    results = []
    file = open(path + '/result/gencode.v44.annotation.gtf')
    line = file.readline()
    while line:
        if '\tgene\t' in line:
            gene = line.split('\t')[8].split(';')[2].split('"')[1]
            if 'ENSG' not in gene:
                results.append(gene)
        line = file.readline()
    return results

def Genion(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'Genion_PacBio.tsv')
    line = PacBio_file.readline()
    while line:
        gene1 = line.split('\t')[1].split('::')[0]
        gene2 = line.split('\t')[1].split('::')[1]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('Genion:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)
    return R

def JAFFAL(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'JAFFAL_PacBio.csv')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'JAFFAL_cDNA.csv')
    line = cDNA_file.readline()
    line = cDNA_file.readline()
    while line:
        if 'Confidence' in line:
            gene1 = line.split(',')[1].split(':')[0]
            gene2 = line.split(',')[1].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    ls = []
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('JAFFAL:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)
    k = []
    for k in new_d[1]:
        if k in new_d[2]:
            ls.append(k)

    return R

def LongGF(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'LongGF_PacBio.log')
    line = PacBio_file.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
    cDNA_file = open(result_dir + 'LongGF_cDNA.log')
    line = cDNA_file.readline()
    while line:
        if 'SumGF' in line:
            gene1 = line.split('\t')[1].split(' ')[0].split(':')[0]
            gene2 = line.split('\t')[1].split(' ')[0].split(':')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l2.append(key)
        line = cDNA_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('LongGF:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)
    return R

def fusionseeker(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'fusionseeker_PacBio.txt')
    line = PacBio_file.readline()
    line = PacBio_file.readline()
    while line:
        gene1 = line.split('\t')[1]
        gene2 = line.split('\t')[2]
        ls = [gene1, gene2]
        ls.sort()
        key = ls[0] + ':' + ls[1]
        l1.append(key)
        line = PacBio_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('fusionseeker:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)
    return R

def GFHunter_RF(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'GFHunter_PacBio.csv')
    line = PacBio_file.readline()
    while line:
        if 'flag = 5' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('GFHunter:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)

    ls = []
    for k in new_d[1]:
        if k in new_d[2]:
            ls.append(k)
    
    return R

def GFHunter(result_dir, result):
    l1 = []
    l2 = []
    l3 = []
    PacBio_file = open(result_dir + 'GFHunter_PacBio.csv')
    line = PacBio_file.readline()
    while line:
        if 'flag = 5' in line or 'flag = 4' in line:
            gene1 = line.split('\t')[0]
            gene2 = line.split('\t')[1]
            ls = [gene1, gene2]
            ls.sort()
            key = ls[0] + ':' + ls[1]
            l1.append(key)
        line = PacBio_file.readline()
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
    new_l1=list(set(l1))
    l1 = new_l1
    new_l2=list(set(l2))
    l2 = new_l2
    new_l3=list(set(l3))
    l3 = new_l3
    l1.sort()
    l2.sort()
    l3.sort()
    data = [l1, l2, l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1 = k.split(':')[0]
            g2 = k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
    print('GFHunter:',end= ' ')
    num = check(new_d)
    R = new_handle(num, data)
    return R

# 新增：计算 ctat-LR-fusion 一致性的函数
# 新增：计算 ctat-LR-fusion 一致性的函数（已增加 num_LR >= 2 过滤）
def ctat(result_dir, result):
    l1 = [] # PB
    l2 = [] # cDNA
    l3 = [] # dRNA
    
    # 读取 PB
    try:
        for ext in ['_PacBio.tsv', '.PacBio.tsv', '_PB.tsv', '.PB.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                line = f.readline()
                line = f.readline()
                while line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        num = int(parts[1]) # 获取支持 reads 数
                        if num >= 2:        # 过滤 >= 2
                            gene1, gene2 = parts[2], parts[5]
                            ls = [gene1, gene2]
                            ls.sort()
                            l1.append(ls[0] + ':' + ls[1])
                    line = f.readline()
                f.close()
                break
    except FileNotFoundError: pass

    # 读取 cDNA
    try:
        for ext in ['_cDNA.tsv', '.cDNA.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                line = f.readline()
                line = f.readline()
                while line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        num = int(parts[1])
                        if num >= 2:
                            gene1, gene2 = parts[2], parts[5]
                            ls = [gene1, gene2]
                            ls.sort()
                            l2.append(ls[0] + ':' + ls[1])
                    line = f.readline()
                f.close()
                break
    except FileNotFoundError: pass

    # 读取 dRNA
    try:
        for ext in ['_dRNA.tsv', '.dRNA.tsv']:
            if os.path.exists(result_dir + 'ctat-LR-fusion' + ext):
                f = open(result_dir + 'ctat-LR-fusion' + ext)
                line = f.readline()
                line = f.readline()
                while line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 8:
                        num = int(parts[1])
                        if num >= 2:
                            gene1, gene2 = parts[2], parts[5]
                            ls = [gene1, gene2]
                            ls.sort()
                            l3.append(ls[0] + ':' + ls[1])
                    line = f.readline()
                f.close()
                break
    except FileNotFoundError: pass

    new_l1, new_l2, new_l3 = list(set(l1)), list(set(l2)), list(set(l3))
    new_l1.sort()
    new_l2.sort()
    new_l3.sort()
    
    data = [new_l1, new_l2, new_l3]
    new_d = [[], [], []]
    for i in range(len(data)):
        for k in data[i]:
            g1, g2 = k.split(':')[0], k.split(':')[1]
            if g1 in result and g2 in result:
                new_d[i].append(k)
                
    print('ctat-LR-fusion:', end=' ')
    num = check(new_d)
    R = new_handle(num, data)
    return R

def new_handle(num, data):
    pcd = num[6]
    cd = num[5] + num[6]
    pd = num[4] + num[6]
    pc = num[3] + num[6]
    pn = len(data[0])
    cn = len(data[1])
    dn = len(data[2])
    
    # 加入防 0 除保护
    r11 = pcd / pn if pn > 0 else 0
    r12 = pcd / cn if cn > 0 else 0
    r13 = pcd / dn if dn > 0 else 0
    R1 = 3 / (1/r11 + 1/r12 + 1/r13) if r11 > 0 and r12 > 0 and r13 > 0 else 0.0
    
    r21 = cd / cn if cn > 0 else 0
    r22 = cd / dn if dn > 0 else 0
    R2 = 2 / (1/r21 + 1/r22) if r21 > 0 and r22 > 0 else 0.0
    
    r31 = pd / pn if pn > 0 else 0
    r32 = pd / dn if dn > 0 else 0
    R4 = 2 / (1/r31 + 1/r32) if r31 > 0 and r32 > 0 else 0.0
    
    r41 = pc / pn if pn > 0 else 0
    r42 = pc / cn if cn > 0 else 0
    R3 = 2 / (1/r41 + 1/r42) if r41 > 0 and r42 > 0 else 0.0

    print(pcd, cd, pd, pc, pn, cn, dn)
    return (R1, R2, R3, R4)


def check(data):
    d = {}
    for i in data[0]:
        if i in d:
            d[i].append(1)
        else:
            d[i] = [1,] 
    for i in data[1]:
        if i in d:
            d[i].append(2)
        else:
            d[i] = [2,]
    for i in data[2]:
        if i in d:
            d[i].append(3)
        else:
            d[i] = [3,]
    num = [0, 0, 0, 0, 0, 0, 0]
    for key,l in d.items():
        if len(l) == 1:
            if l[0] == 1:
                num[0] += 1
            elif l[0] == 2:
                num[1] += 1
            else:
                num[2] += 1
        elif len(l) == 2:
            if 3 not in l:
                num[3] += 1
            elif 2 not in l:
                num[4] += 1
            else:
                num[5] += 1
        else:
            num[6] += 1
    #print(num)
    return num

if __name__ == "__main__":
    r = read_result()
    dir = path + '/result/MCF-7/'
    R = []
    R.append(GFHunter_RF(dir, r))
    R.append(GFHunter(dir, r))
    R.append(LongGF(dir, r))
    R.append(JAFFAL(dir, r))
    R.append(fusionseeker(dir, r))
    R.append(Genion(dir, r))
    R.append(ctat(dir, r)) # 调用 ctat 获取结果

    with open(outputdir + 'ratioupset.csv', 'w') as f:
        f.write('tools\tcDNA+dRNA+PB\tcDNA+dRNA\tcDNA+PB\tdRNA+PB\n')
        line1 = 'GFHunter RF\t' + str(R[0][0]) + '\t' +  str(R[0][1]) + '\t' +  str(R[0][2]) + '\t' +  str(R[0][3]) + '\n'
        line2 = 'GFHunter RF+SF\t' + str(R[1][0]) + '\t' +  str(R[1][1]) + '\t' +  str(R[1][2]) + '\t' +  str(R[1][3]) + '\n'
        line3 = 'LongGF\t' + str(R[2][0]) + '\t' +  str(R[2][1]) + '\t' +  str(R[2][2]) + '\t' +  str(R[2][3]) + '\n'
        line4 = 'JAFFAL HC+LC\t' + str(R[3][0]) + '\t' +  str(R[3][1]) + '\t' +  str(R[3][2]) + '\t' +  str(R[3][3]) + '\n'
        line5 = 'fusionseeker\t' + str(R[4][0]) + '\t' +  str(R[4][1]) + '\t' +  str(R[4][2]) + '\t' +  str(R[4][3]) + '\n'
        line6 = 'Genion\t' + str(R[5][0]) + '\t' +  str(R[5][1]) + '\t' +  str(R[5][2]) + '\t' +  str(R[5][3]) + '\n'
        line7 = 'ctat-LR-fusion\t' + str(R[6][0]) + '\t' +  str(R[6][1]) + '\t' +  str(R[6][2]) + '\t' +  str(R[6][3]) + '\n'
        
        f.write(line1)
        f.write(line2)
        f.write(line3)
        f.write(line4)
        f.write(line5)
        f.write(line6)
        f.write(line7) # 写入 ctat 数据行