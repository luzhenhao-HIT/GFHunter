import time
import numpy as np

K = 200

#annotation = '../middlefile/rebuild_annotation.txt'
sam_file = 'middlefile/step4.RNA_alignment.only_primary.sam'
result_in = 'middlefile/step3.genefusion.csv'
result_out = 'step5.genefusion.csv'

'''
class: Gene
'''
class Gene():
    def __init__(self, id, name, start, end, chr, flag):
        self.id = id
        self.name = name
        self.start = start
        self.end = end
        self.chr = chr
        self.flag = flag
        pass
'''
class: Fusion
'''
class Fusion():
    def __init__(self, gene1, gene2, num, bp1, bp2, code):
        self.gene1 = gene1
        self.gene2 = gene2
        self.num = num
        self.bp1 = bp1
        self.bp2 = bp2
        self.code = code
        self.read_id = self.gene1.id + '_' + self.gene2.id + '_' + code
        self.flag = False
        self.samflag = False
        pass
    def add_samlines(self, samlines):
        '''if self.gene1.id == 'ENSG00000170421.13':
            print(samlines)'''
        self.samlines = samlines
        if samlines != []:
            self.samflag = True
        self.sort_samlines()
        pass
    def sort_samlines(self):
        read_to_samlines = {}
        true_read = []
        self.bp1s = []
        self.bp2s = []
        for samline in self.samlines:
            if samline.read in read_to_samlines:
                read_to_samlines[samline.read].append(samline)
            else:
                read_to_samlines[samline.read] = [samline,]
        for read, samlines in read_to_samlines.items():
            samlines1 = []
            samlines2 = []
            flag = False
            for samline in samlines:
                if samline.end - self.gene1.start > 0 and self.gene1.end - samline.start > 0:
                    bp1, f = self.check_breakpoints1(samline)
                    '''if self.gene1.id == 'ENSG00000170421.13':
                        print(bp1, f)'''
                    if  f == True:
                        samlines1.append(samline)
                        self.bp1s.append(bp1)
                elif samline.end - self.gene2.start > 0 and self.gene2.end - samline.start > 0:
                    bp2, f = self.check_breakpoints2(samline)
                    '''if self.gene1.id == 'ENSG00000170421.13':
                        print(bp2,f)'''
                    if  f == True:
                        samlines2.append(samline)
                        self.bp2s.append(bp2)
            if samlines1 != [] and samlines2 != []:
                for s1 in samlines1:
                    for s2 in samlines2:
                        if s1.position1 < s2.position1:
                            flag = True
                            self.flag = True
                            break
                    if flag == True:
                        true_read.append(read)
                        break
        n = 0
        for read in true_read:
            n += int(read.split('_')[3].split('\n')[0])
        if self.flag == True:
            self.newnum = n
            self.bp1 = int(np.mean(self.bp1s))
            self.bp2 = int(np.mean(self.bp2s))
        '''if self.gene1.id == 'ENSG00000170421.13':
            print(self.bp1, self.bp2, self.flag, self.samflag)'''
        pass
    def check_breakpoints1(self, samline):
        flag = False
        bp = 0
        if samline.flag == 0 or samline.flag == 256 or samline.flag == 2048:
            if np.abs(samline.end - self.bp1) <= K:
                flag = True
                bp = samline.end
        elif samline.flag == 16 or samline.flag == 272 or samline.flag == 2064:
            if np.abs(samline.start - self.bp1) <= K:
                flag = True
                bp = samline.start
        return(bp, flag)
        pass
    def check_breakpoints2(self, samline):
        flag = False
        bp = 0
        if samline.flag == 0 or samline.flag == 256 or samline.flag == 2048:
            if np.abs(samline.start - self.bp2) <= K:
                flag = True
                bp = samline.start
        elif samline.flag == 16 or samline.flag == 272 or samline.flag == 2064:
            if np.abs(samline.end - self.bp2) <= K:
                flag = True
                bp = samline.end
        return (bp, flag)
        pass
'''
class: samline
'''
class Samline():
    def __init__(self, read, chr, flag, start, cigar, num):
        self.read = read
        self.chr = chr
        self.start = start
        self.cigar = cigar
        self.flag = flag
        self.num = num
        self.compute_end()
        self.compute_position()
        pass
    def compute_end(self):
        num = ''
        length = 0
        for i in self.cigar:
            if i.isdigit():
                num = num + i
            elif i == 'I' or i == 'S' or i == 'H':
                num = ''
            else:
                length = length + int(num)
                num = ''
        self.length = length
        self.end = self.start + length - 1
        '''if 'ENSG00000170421.13' in self.read:
            print(self.start, self.end)'''
        pass
    def compute_position(self):
        num = ''
        length = 0
        for i in self.cigar:
            if i.isdigit():
                num = num + i
            elif i == 'D'  or i == 'N':
                num = ''
            else:
                length = length + int(num)
                num = ''
        self.ownlength = length
        if 'S' in self.cigar or 'H' in self.cigar:
            if 'S' in self.cigar:
                key = 'S'
            elif 'H' in self.cigar:
                key = 'H'
            num = len(self.cigar.split(key))
            if num == 2:
                if self.cigar.split(key)[1] == '':
                    self.position1 = 0
                    num = ''
                    for i in self.cigar.split(key)[0][::-1]:
                        if i.isdigit():
                            num += i
                        else:
                            break
                    num = int(num[::-1])
                    self.position2 = length - num
                else:
                    self.position1 = int(self.cigar.split(key)[0])
                    self.position2 = length
            else:
                self.position1 = int(self.cigar.split(key)[0])
                num = ''
                for i in self.cigar.split(key)[1][::-1]:
                    if i.isdigit():
                        num += i
                    else:
                        break
                num = int(num[::-1])
                self.position2 = length - num
        else:
            self.position1 = 0
            self.position2 = length
        if self.flag == 16 or self.flag == 2064:
            position1 = self.ownlength - self.position1 + 1
            position2 = self.ownlength - self.position2 + 1
            self.position1 = position2
            self.position2 = position1
        '''if 'ENSG00000170421.13' in self.read:
            print(self.position1, self.position2)'''
'''
fuction: read_result
'''
def read_result(result_in):
    fusions = {}
    file = open(result_in)
    line = file.readline()
    line = file.readline()
    while line:
        n1 = line.split('\t')[0]
        n2 = line.split('\t')[1]
        id1 = line.split('\t')[2]
        id2 = line.split('\t')[3]
        chr1 = line.split('\t')[4].split(':')[0]
        flag1 = line.split('\t')[4].split(':')[1]
        chr2 = line.split('\t')[5].split(':')[0]
        flag2 = line.split('\t')[5].split(':')[1]
        start1 = int(line.split('\t')[6].split(':')[0])
        end1 = int(line.split('\t')[6].split(':')[1])
        start2 = int(line.split('\t')[7].split(':')[0])
        end2 = int(line.split('\t')[7].split(':')[1])
        bp1 = int(line.split('\t')[8])
        bp2 = int(line.split('\t')[9])
        num = int(line.split('\t')[10])
        code = line.split('\t')[11].split('\n')[0]

        gene1 = Gene(id1, n1, start1, end1, chr1, flag1)
        gene2 = Gene(id2, n2, start2, end2, chr2, flag2)
        fusion = Fusion(gene1, gene2, num, bp1, bp2, code)
        if id1 + '_' + id2 in fusions:
            fusions[id1 + '_' + id2].append(fusion)
        else:
            fusions[id1 + '_' + id2] = [fusion,]
        line = file.readline()
    return fusions
'''
fuction: read_sam_file
'''
def read_sam_file(sam_file):
    read_to_samlines_with_SA = {} #{key = read : value = [samline,...]}
    read_to_samline_without_SA = {} #{key = read : value = samline}
    samfile = open(sam_file)
    line = samfile.readline()
    while line:
        if '@' not in line:
            flag = int(line.split('\t')[1])
            if flag != 4:
                read = line.split('\t')[0]
                id1 = read.split('_')[0]
                id2 = read.split('_')[1]
                read_id = id1 + '_' + id2
                chr = line.split('\t')[2]
                start = int(line.split('\t')[3])
                cigar = line.split('\t')[5]
                num = read.split('_')[3].split('\n')[0]
                samline = Samline(read, chr, flag, start, cigar, num)
                if read_id in read_to_samlines_with_SA:
                    read_to_samlines_with_SA[read_id].append(samline)
                else:
                    read_to_samlines_with_SA[read_id] = [samline,]
                if 'SA' in line:
                    line = samfile.readline()
                    while line.split('\t')[0] == read:
                        read = line.split('\t')[0]
                        chr = line.split('\t')[2]
                        start = int(line.split('\t')[3])
                        cigar = line.split('\t')[5]
                        flag = int(line.split('\t')[1])
                        num = read.split('_')[3].split('\n')[0]
                        samline = Samline(read, chr, flag, start, cigar, num)
                        read_to_samlines_with_SA[read_id].append(samline)
                        line = samfile.readline()
                else:
                    line = samfile.readline()
            else:
                line = samfile.readline()
        else:
            line = samfile.readline()
    return read_to_samlines_with_SA
'''
fuction: outputline
'''
def outputline(fusion):
    l1 = fusion.gene1.name
    l2 = fusion.gene2.name
    l3 = fusion.gene1.id
    l4 = fusion.gene2.id
    l5 = fusion.gene1.chr + ':' + str(fusion.bp1) + '; ' + fusion.gene2.chr + ':' + str(fusion.bp2)
    l6 = fusion.gene1.chr + fusion.gene1.flag + ': ' + str(fusion.gene1.start) + ', ' + str(fusion.gene1.end)
    l7 = fusion.gene2.chr + fusion.gene2.flag + ': ' + str(fusion.gene2.start) + ', ' + str(fusion.gene2.end)
    if fusion.samflag == False:
        flag = 3
        l8 = str(fusion.num)
        l9 = 'flag = 3: Filtering'
    elif fusion.samflag == True and fusion.flag == False:
        flag = 4
        l8 = str(fusion.num)
        l9 = 'flag = 4: Pseudo fusion gene'
    elif fusion.flag == True:
        flag = 5
        l8 = str(fusion.newnum)
        l9 = 'flag = 5: Fusion gene'
    line = l1 + '\t' + l2 + '\t' + l3 + '\t' + l4 + '\t' + l5 + '\t' + l6 + '\t' + l7 + '\t' + l8 + '\t' + l9 + '\n'
    return line
        
'''
fuction: filtering_and_output
'''
def filtering_and_output(read_to_samline, fusions, result_out):
    with open(result_out, 'w') as file:
        line = 'gene1 name\tgene2 name\tgene1 id\tgene2 id\tbreak points\tgene1 position\tgene2 position\tsupport num\tflag\n'
        file.write(line)
        for read_id, fusions in fusions.items():
            samlines = []
            if read_id in read_to_samline:
                samlines = read_to_samline[read_id]
            for fusion in fusions:
                fusion.add_samlines(samlines)
                line = outputline(fusion)
                file.write(line)
def main(result_out, middlefile):
    sam_file = middlefile + 'step2.RNA_alignment.only_primary.sam'
    result_in = middlefile + 'step1.genefusion.csv'
    fusions = read_result(result_in)
    read_to_samline = read_sam_file(sam_file)
    filtering_and_output(read_to_samline, fusions, result_out)

if __name__ == "__main__":
    fusions = read_result(result_in)
    read_to_samline = read_sam_file(sam_file)
    filtering_and_output(read_to_samline, fusions, result_out)