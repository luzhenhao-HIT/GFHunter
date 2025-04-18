import time, sys, os
import shlex
import numpy as np
from subprocess import check_output, PIPE, Popen

K = 200
least_reads = 2
POA = ''

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
    def __init__(self, gene1, gene2, num, bp1, bp2, code, reads, only_id):
        self.gene1 = gene1
        self.gene2 = gene2
        self.num = num
        self.bp1 = bp1
        self.bp2 = bp2
        self.code = code
        self.reads = reads
        self.read_id = self.gene1.id + '_' + self.gene2.id + '_' + code
        self.only_id = only_id
        self.flag = 2
        pass
    def add_result(self, result):
        self.result = result
        self.flag = 3
'''
class: Result
'''
class Result():
    def __init__(self, only_id, sequence, samline):
        self.only_id = only_id
        self.transcript = sequence
        self.samlines = [samline,]
        self.num = samline.num
        self.flag = 3
        self.reads = ''
        pass
    def add_samline(self, samline):
        self.samlines.append(samline)
        pass
    def add_fusion(self, fusion):
        self.fusion = fusion
        self.reads = fusion.reads
        pass
    def add_bp(self, bp1s, bp2s):
        self.bp1 = int(np.mean(bp1s))
        self.bp2 = int(np.mean(bp2s))
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
        if self.flag == 16 or self.flag == 2064 or self.flag == 272:
            position1 = self.ownlength - self.position1 + 1
            position2 = self.ownlength - self.position2 + 1
            self.position1 = position2
            self.position2 = position1
        pass

'''
fuction: Antisense
input: string
output: antistring
'''
def Antisense(string):
    antistring = ''
    for i in string:
        if i == 'A':
            antistring += 'T'
        elif i == 'T':
            antistring += 'A'
        elif i == 'C':
            antistring += 'G'
        elif i == 'G':
            antistring += 'C'
    antistring = antistring[::-1]
    return antistring
'''
fuction: subprocess_popen
'''
def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=sys.stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)
'''
fuction: read_result
'''
def read_result(result_file):
    fusions = {}
    for line in result_file.split('\n')[:-1]:
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
        code = line.split('\t')[11]
        reads = line.split('\t')[12].split('\n')[0]
        only_id = id1 + '_' + id2 + '_' + code + '_' + str(num)
        gene1 = Gene(id1, n1, start1, end1, chr1, flag1)
        gene2 = Gene(id2, n2, start2, end2, chr2, flag2)
        fusion = Fusion(gene1, gene2, num, bp1, bp2, code, reads, only_id)
        if id1 + '_' + id2 in fusions:
            fusions[id1 + '_' + id2].append(fusion)
        else:
            fusions[id1 + '_' + id2] = [fusion,]
    return fusions
def read_result0():
    file = open('./middlefile/step1.genefusion.csv')
    line = file.readline()
    linelist = file.readlines()
    result = ''
    for line in linelist:
        result += line
    return result

'''
fuction: alignment
'''
def alignment(POA_file, index_dir, middlefile, thread, if_output):
    results_of_each_fusion = {} #{key = key_id : value = [result,...]}
    ref_index = index_dir + 'ref_index.mmi'
    minimap2_view_process = subprocess_popen(shlex.split('minimap2 -ax splice -uf -k14  -t {} {} {}'.format(str(thread), ref_index, POA_file)))
    with open(middlefile + 'step2.RNA_alignment.sam', 'w') as samfile:
        read0 = ''
        for line in minimap2_view_process.stdout:
            if line != '' and line != None:
                if if_output == True:
                    samfile.write(line)
                if '@' not in line:
                    flag = int(line.split('\t')[1])
                    if flag != 4:
                        read = line.split('\t')[0]
                        if read != read0:
                            read0 = read
                            id1 = read.split('_')[0]
                            id2 = read.split('_')[1]
                            num = int(read.split('_')[3])
                            read_id = id1 + '_' + id2
                            chr = line.split('\t')[2]
                            start = int(line.split('\t')[3])
                            cigar = line.split('\t')[5]
                            if flag == 0:
                                sequence = line.split('\t')[9]
                            elif flag == 16:
                                sequence = Antisense(line.split('\t')[9])
                            samline = Samline(read, chr, flag, start, cigar, num)
                            result = Result(read, sequence, samline)
                            if read_id in results_of_each_fusion:
                                results_of_each_fusion[read_id].append(result)
                            else:
                                results_of_each_fusion[read_id] = [result,]
                        else:
                            id1 = read.split('_')[0]
                            id2 = read.split('_')[1]
                            num = int(read.split('_')[3])
                            read_id = id1 + '_' + id2
                            chr = line.split('\t')[2]
                            start = int(line.split('\t')[3])
                            cigar = line.split('\t')[5]
                            samline = Samline(read, chr, flag, start, cigar, num)
                            result.add_samline(samline)
    return results_of_each_fusion

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
    l10 = fusion.reads
    if fusion.flag == 2:
        flag = 3
        l8 = str(fusion.num)
        l9 = 'flag = 3: Potential Fusion'
    if fusion.flag == 3:
        if fusion.result.flag == 3:
            flag = 3
            l8 = str(fusion.num)
            l9 = 'flag = 3: Potential Fusion'
        elif fusion.result.flag == 4:
            flag = 4
            l8 = str(fusion.num)
            l9 = 'flag = 4: Suspected Fusion'
        elif fusion.result.flag == 5:
            flag = 5
            l5 = fusion.gene1.chr + ':' + str(fusion.result.bp1) + '; ' + fusion.gene2.chr + ':' + str(fusion.result.bp2)
            l8 = str(fusion.num)
            l9 = 'flag = 5: Reliable Fusion'
    line = l1 + '\t' + l2 + '\t' + l3 + '\t' + l4 + '\t' + l5 + '\t' + l6 + '\t' + l7 + '\t' + l8 + '\t' + l9 + '\t' + l10 + '\n'
    return flag, line
        
'''
fuction: filtering_and_output
'''
def filtering_and_output(results_of_each_fusion, fusions_dict, result_out):
    with open(result_out, 'w') as file:
        line = 'gene1 name\tgene2 name\tgene1 id\tgene2 id\tbreak points\tgene1 position\tgene2 position\tsupport num\trank class\tsupperting reads information\n'
        file.write(line)
        for key, fusions in fusions_dict.items():
            if key in results_of_each_fusion:
                gene1 = fusions[0].gene1
                gene2 = fusions[0].gene2
                results = results_of_each_fusion[key]
                for result in results:
                    sam1 = []
                    sam2 = []
                    for samline in result.samlines:
                        if samline.end - gene1.start > 0 and gene1.end - samline.start > 0 and samline.chr == gene1.chr:
                            sam1.append(samline)
                        elif samline.end - gene2.start > 0 and gene2.end - samline.start > 0 and samline.chr == gene2.chr:
                            sam2.append(samline)
                    if sam1 != [] and sam2 != []:
                        result.flag = 4
                        for fusion in fusions:
                            if result.only_id == fusion.only_id:
                                result.add_fusion(fusion)
                                fusion.add_result(result)
                                break
                        for fusion in fusions:
                            bp1s = []
                            bp2s = []
                            for sam in sam1:
                                bp1, f1 = check_breakpoints1(fusion.bp1, sam)
                                if f1 == True:
                                    bp1s.append(bp1)
                            for sam in sam2:
                                bp2, f2 = check_breakpoints2(fusion.bp2, sam)
                                if f2 == True:
                                    bp2s.append(bp2)
                            if bp1s != [] and bp2s != []:
                                result.add_bp(bp1s, bp2s)
                                if result.only_id == fusion.only_id:
                                    result.flag = 5
                                else:
                                    fusion.num += result.num
                                    fusion.reads += result.reads
        for key, fusions in fusions_dict.items():
            for fusion in fusions:
                if (fusion.flag == 2 and fusion.num >= least_reads) or (fusion.flag == 3 and fusion.result.num >= least_reads):
                    flag, line = outputline(fusion)
                    file.write(line)
        '''lines_5 = []
        lines_4 = []
        lines_3 = []
        for key, fusions in fusions_dict.items():
            for fusion in fusions:
                if (fusion.flag == 2 and fusion.num >= least_reads) or (fusion.flag == 3 and fusion.result.num >= least_reads):
                    flag, line = outputline(fusion)
                    if flag == 5:
                        lines_5.append(line)
                    elif flag == 4:
                        lines_4.append(line)
                    elif flag == 3:
                        lines_3.append(line)
        for line in lines_5:            
            file.write(line)
        for line in lines_4:            
            file.write(line)
        for line in lines_3:            
            file.write(line)'''

def check_breakpoints1(bp1, samline):
    flag = False
    bp = 0
    if samline.flag == 0 or samline.flag == 256 or samline.flag == 2048:
        if np.abs(samline.end - bp1) <= K:
            flag = True
            bp = samline.end
    elif samline.flag == 16 or samline.flag == 272 or samline.flag == 2064:
        if np.abs(samline.start - bp1) <= K:
            flag = True
            bp = samline.start
    return (bp, flag)

def check_breakpoints2(bp2, samline):
    flag = False
    bp = 0
    if samline.flag == 0 or samline.flag == 256 or samline.flag == 2048:
        if np.abs(samline.start - bp2) <= K:
            flag = True
            bp = samline.start
    elif samline.flag == 16 or samline.flag == 272 or samline.flag == 2064:
        if np.abs(samline.end - bp2) <= K:
            flag = True
            bp = samline.end
    return (bp, flag)


def main(result_out, POA_file, result_file, index_dir, middlefile, thread, if_output, lreads, bpsbp):
    global least_reads, K
    K = bpsbp
    least_reads = lreads
    fusions = read_result(result_file)
    results_of_each_fusion = alignment(POA_file, index_dir, middlefile, thread, if_output)
    filtering_and_output(results_of_each_fusion, fusions, result_out)

if __name__ == "__main__":
    result_out = './result.csv'
    POA_file = './middlefile/step1.POA_result.fasta'
    result_file = read_result0()
    index_dir = '../GFHunter/index/'
    middlefile = './middlefile/'
    thread = 4
    if_output = True
    main(result_out, POA_file, result_file, index_dir, middlefile, thread, if_output, 2)
    