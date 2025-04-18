import time
import os, sys
import shlex
from subprocess import check_output, PIPE, Popen
from intervaltree import IntervalTree
import numpy as np
from scipy.sparse import csr_matrix
import pyabpoa as pa

min_length = 50
exon_boundary = 30
overlap_present = 0.5
clustering_num = 200
sup_read_num = 2
type_alignment = 'ont'
sam_file = 'middlefile/step1.DNA_like_alignment.only_primary.sam'
output_file_dir = 'middlefile/step1.gene_groups/'
tree_dir = 'middlefile/step1.clustering_tree/'
out_put = 'middlefile/step1.POA_result.fasta'
result_out = 'middlefile/step1.genefusion.csv'
output_sam = 'middlefile/step1.DNA_like_alignment.sam'

if_output = True
thread = 32
min_abpoa = 100
'''
global variable
'''
index = {} #{key = chr : value = IntervalTree_of_transcript_in_the_chr}
reads = [] #[read1, read2, ...]  readn = [read_tuple, ...]
'''
class: Gene
'''
class Gene():
    def __init__(self, id, name, refstart, refend, chr, flag, gene_type):
        self.id = id
        self.name = name
        self.rstart = refstart
        self.rend = refend
        self.chr = chr
        self.flag = flag
        self.type = gene_type
        self.transcripts = []
        pass
    def add_transcript(self, transcript):
        self.transcripts.append(transcript)
        pass
'''
class: Transcript
'''
class Transcript():
    def __init__(self, id, name, refstart, refend, start, end, chr, flag, trans_type, gene):
        self.id = id
        self.name = name
        self.rstart = refstart
        self.rend = refend
        self.start = start
        self.end = end
        self.chr = chr
        self.flag = flag
        self.trans_type = trans_type
        self.gene = gene
        self.exons = []
        pass
    def add_exon(self, exon):
        self.exons.append(exon)
        pass
'''
class: Exon
'''
class Exon():
    def __init__(self, id, refstart, refend, start, end, transcript):
        self.id = id
        self.rstart = refstart
        self.rend = refend
        self.start = start
        self.end = end
        self.transcript = transcript
        self.length = self.end - self.start + 1
        pass
'''
class Samline
'''
class Samline():
    def __init__(self, id, start, cigar, flag, chr, fasta):
        self.id = id
        self.start = start
        self.cigar = cigar
        self.flag = flag
        self.chr = chr
        self.transcripts = []
        self.overlap = 0
        self.transcript = 0
        self.sequence = fasta
        if flag == 0 or flag == 16:
            self.major = True
        else:
            self.major = False
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
            elif i == 'D':
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
            if 'S' in self.cigar:
                self.sequence = self.sequence[self.position1:self.position2]
        else:
            self.position1 = 0
            self.position2 = length
        if self.flag == 16 or self.flag == 2064:
            position1 = self.ownlength - self.position1 + 1
            position2 = self.ownlength - self.position2 + 1
            self.position1 = position2
            self.position2 = position1
            self.sequence = Antisense(self.sequence)
    def add_num(self, num):
        self.num = num
    def add_transcript(self, transcript, start, end):
        self.transcripts.append(transcript)
        overlap1 = end - self.start
        overlap2 = self.end - start
        length1 = self.length
        length2 = end - start +1
        real_overlap = min(overlap1, overlap2, length1, length2)
        if self.overlap <= real_overlap:
            self.overlap = real_overlap
            self.transcript = transcript
            self.tstart = start
            self.tend = end
        pass
    def add_read(self, read):
        self.read = read
'''
class: Read
'''
class Read():
    def __init__(self, id, start, cigar, flag, chr, fasta):
        self.id = id
        self.num = 1
        self.samlines = []
        self.sequence = fasta
        if flag == 16:
            self.sequence = Antisense(self.sequence)
        samline = Samline(id, start, cigar, flag, chr, fasta)
        samline.compute_end()
        samline.compute_position()
        samline.add_read(self)
        length = samline.position2 - samline.position1 + 1
        if length > min_length:
            self.samlines.append(samline)
        pass
    def add_samline(self, id, start, cigar, flag, chr, fasta):
        self.num += 1
        samline = Samline(id, start, cigar, flag, chr, fasta)
        samline.compute_end()
        samline.compute_position()
        samline.add_read(self)
        length = samline.position2 - samline.position1 + 1
        if length > min_length:
            self.samlines.append(samline)
        pass
    def num_samline(self):
        samlines = self.samlines
        n = len(samlines)
        for i in range(n - 1):
            for j in range(i + 1, n):
                if samlines[i].position1 > samlines[j].position1:
                    samline = samlines[j]
                    samlines[j] = samlines[i]
                    samlines[i] = samline
        for i in range(n):
            samlines[i].add_num(i)
        pass
    def transcripts_and_genes(self):
        global exon_boundary
        self.transcripts = []
        self.genes = {}
        samlines = []
        self.read_through = []
        poslist = []
        for samline in self.samlines:
            if samline.transcript != 0:                
                if samline.start > samline.transcript.start + exon_boundary or samline.end < samline.transcript.end - exon_boundary: #?????
                    samlines.append(samline)
                else:
                    self.read_through.append(samline)
        if samlines == []:
            return False
        self.samlines = samlines
        for samline in self.samlines:
            if samline.transcript not in self.transcripts:
                self.transcripts.append(samline.transcript)
                if samline.transcript.gene not in self.genes:
                    self.genes[samline.transcript.gene] = [samline.transcript,]
                else:
                    self.genes[samline.transcript.gene].append(samline.transcript)
        return True
    def group_samlines(self, gene1, gene2):
        length1 = 0
        length2 = 0
        for samline in self.samlines:
            if samline.transcript.gene.id == gene1.id:
                if samline.length > length1:
                    length1 = samline.length
                    samline1 = samline
            elif samline.transcript.gene.id == gene2.id:
                if samline.length > length2:
                    length2 = samline.length
                    samline2 = samline
        return (samline1, samline2)
'''
class: Genes_group
'''
class Genes_group():
    def __init__(self, gene1, gene2, reads):
        self.gene1 = gene1
        self.gene2 = gene2
        self.reads = reads
        pass
    def clustering(self):
        global sup_read_num
        self.trees = {}
        self.leaves = {}
        self.flag = False
        n = 0
        for read in self.reads:
            samline1, samline2 = read.group_samlines(self.gene1, self.gene2)
            leaf = Tree_leaf(read, n, samline1, samline2)
            if leaf.flag != 4:
                if leaf.flag in self.leaves:
                    self.leaves[leaf.flag].append(leaf)
                else:
                    self.leaves[leaf.flag] = [leaf,]
                n += 1
        for flag, leaves in self.leaves.items():
            if len(leaves) >= sup_read_num:
                self.flag = True
                self.trees[flag] = Cluster_tree(leaves, self.gene1, self.gene2, flag)
                self.trees[flag].clustering() 
        pass
'''
class: Cluster_tree
'''
class Tree_leaf():
    def __init__(self, read, num, samline1, samline2):
        self.read = read
        self.num = num
        self.samline1 = samline1
        self.samline2 = samline2
        self.bp_position()
        pass
    def add_address(self, address):
        self.address = address
        pass
    def front_bp_pos(self, samline):
        global exon_boundary
        for exon in samline.transcript.exons[:-1]:
            if exon.end - exon_boundary <= samline.end and exon.end + exon_boundary >= samline.end:
                return samline.end - exon.end + exon.rend
        return False
    def behind_bp_pos(self, samline):
        global exon_boundary
        for exon in samline.transcript.exons[1:]:
            if exon.start - exon_boundary <= samline.start and exon.start + exon_boundary >= samline.start:
                return samline.start - exon.start + exon.rstart
        return False
    def bp_position(self):
        if self.samline1.num < self.samline2.num:
            #gene1 -> gene2
            front_samline = self.samline1
            behind_samline = self.samline2
            self.anti = False
        else:
            #gene2 -> gene1
            front_samline = self.samline2
            behind_samline = self.samline1
            self.anti = True
        self.fgene = front_samline.transcript.gene
        self.bgene = behind_samline.transcript.gene
        self.fsam = front_samline
        self.bsam = behind_samline
        if front_samline.flag == 0 or front_samline.flag == 2048:
            if behind_samline.flag == 0 or behind_samline.flag == 2048:
                #++
                self.bp1 = self.front_bp_pos(front_samline)
                self.bp2 = self.behind_bp_pos(behind_samline)
                self.flag = 0
                pass
            else:
                #+-
                self.bp1 = self.front_bp_pos(front_samline)
                self.bp2 = self.front_bp_pos(behind_samline)
                self.flag = 1
                pass
        else:
            if behind_samline.flag == 0 or behind_samline.flag == 2048:
                #-+
                self.bp1 = self.behind_bp_pos(front_samline)
                self.bp2 = self.behind_bp_pos(behind_samline)
                self.flag = 2
                pass
            else:
                #--
                self.bp1 = self.behind_bp_pos(front_samline)
                self.bp2 = self.front_bp_pos(behind_samline)
                self.flag = 3
                pass
        if self.samline1.num > self.samline2.num:
            if self.flag == 0 or self.flag == 3:
                self.flag = 3 - self.flag
        if self.bp1 == False or self.bp2 == False:
            self.flag = 4
class Tree_node():
    def __init__(self, num, height, leaves):
        self.num = num
        self.height = height
        self.leaves = leaves
        self.father = None
        self.leftson = None
        self.rightson = None
        self.calculate_address()
        pass
    def add_group(self, group, group_flag):
        self.group = group
        self.gf = group_flag
        pass
    def calculate_address(self):
        address1 = 0
        address2 = 0
        for leaf in self.leaves:
            address1 += leaf.address[0]
            address2 += leaf.address[1]
        address1 = int(address1/len(self.leaves))
        address2 = int(address2/len(self.leaves))
        self.address = (address1, address2)
    def add_father(self, father):
        self.father = father
    def add_sons(self, son1, son2):
        self.leftson = son1
        self.rightson = son2
class Cluster_tree():
    def __init__(self, leaves, gene1, gene2, flag):
        self.leaflist = leaves
        self.flag = flag
        address1list1 = []
        address1list2 = []
        self.root = None
        self.nodes = {0 : []}
        for leaf in self.leaflist:
            if leaf.fgene == gene1:
                address1list1.append(leaf.bp1)
                address1list2.append(leaf.bp2)
            else:
                address1list1.append(leaf.bp2)
                address1list2.append(leaf.bp1)
        mean1 = np.mean(address1list1)
        mean2 = np.mean(address1list2)
        for i in range(len(self.leaflist)):
            leaf = self.leaflist[i]
            leaf.num = i
            leaf.add_address((address1list1[i] - mean1, address1list2[i] - mean2))
        num = 0
        for leaf in self.leaflist:
            node = Tree_node(num, 0, [leaf,])
            num += 1
            self.nodes[0].append(node)
        self.num = num - 1
        self.mean1 = int(mean1)
        self.mean2 = int(mean2)
        pass
    def clustering(self):
        n = 0
        while len(self.nodes[n]) > 1:
            self.nodes[n + 1] = self.find_min_distance(self.nodes[n])
            n += 1
        self.root = self.nodes[n][0]
        pass
    def find_min_distance(self, nodes):
        n = len(nodes)
        distance = self.calculate_distance(nodes[0].address, nodes[1].address)
        lnode = nodes[0]
        rnode = nodes[1]
        for i in range(n - 1):
            for j in range(i + 1, n):
                if distance > self.calculate_distance(nodes[i].address, nodes[j].address):
                    distance = self.calculate_distance(nodes[i].address, nodes[j].address)
                    lnode = nodes[i]
                    rnode = nodes[j]
        self.num += 1
        leaves = lnode.leaves + rnode.leaves
        fnode = Tree_node(self.num, distance, leaves)
        fnode.add_sons(lnode, rnode)
        lnode.add_father(fnode)
        rnode.add_father(fnode)
        newnodes = [fnode,]
        for node in nodes:
            if node != rnode and node != lnode:
                newnodes.append(node)
        return newnodes
    def calculate_distance(self, a, b):
        distance = abs(a[0] - b[0]) + abs(a[1] - b[1])
        return distance
    def print_tree(self, file):
        with open(file, 'w') as result:
            for i in range(len(self.nodes.keys())):
                n = len(self.nodes.keys()) - i - 1
                line0 = 'layer' + str(n) + ':\n'
                line1 = ''
                for node in self.nodes[n]:
                    line1 += 'node num: ' + str(node.num) + '\tdistance between sons: ' + str(node.height) + '\n' + 'leaf num: '
                    for leaf in node.leaves:
                        line1 += str(leaf.num) + ' '
                    line1 += '\naddress: ' + str(node.address[0] + self.mean1) + ', ' + str(node.address[1] + self.mean2) + '\n'
                result.write(line0)
                result.write(line1 + '\n')
        pass
    def classify(self, k):
        global sup_read_num
        for i in range(len(self.nodes.keys())):
            n = len(self.nodes.keys()) - i - 1
            if_classify = True
            for node in self.nodes[n]:
                if node.height > k:
                    if_classify = False
            if if_classify == True :
                min = sup_read_num
                nodelist = []
                for node in self.nodes[n]:
                    if len(node.leaves) >= min:
                        nodelist.append(node)
                return nodelist
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
step 1.1
fuction: read_annotation()
input: annotation
output: NULL (global index)
Read rebuilded annotation to make index for samlines to gene.
'''
def read_annotation(annotation):
    global index
    m = 0
    n = 0
    chr_to_genes = {}
    gtf = open(annotation)
    line = gtf.readline()
    line = gtf.readline()
    while line:
        if line.split('\t')[1] == 'gene':
            m += 1
            id = line.split('\t')[5]
            name = line.split('\t')[6]
            refstart = int(line.split('\t')[2].split(':')[0])
            refend = int(line.split('\t')[2].split(':')[1])
            chr = line.split('\t')[0]
            flag = line.split('\t')[4]
            gene_type = line.split('\t')[7].split('\n')[0]
            gene = Gene(id, name, refstart, refend, chr, flag, gene_type)
            if chr in chr_to_genes:
                chr_to_genes[chr].append(gene)
            else:
                chr_to_genes[chr] = [gene,]
            pass
        elif line.split('\t')[1] ==  'transcript':
            n += 1
            id = line.split('\t')[5]
            name = line.split('\t')[6]
            refstart = int(line.split('\t')[2].split(':')[0])
            refend = int(line.split('\t')[2].split(':')[1])
            start = int(line.split('\t')[3].split(':')[0])
            end = int(line.split('\t')[3].split(':')[1])
            chr = line.split('\t')[0]
            flag = line.split('\t')[4]
            trans_type = line.split('\t')[7].split('\n')[0]
            transcript = Transcript(id, name, refstart, refend, start, end, chr, flag, trans_type, gene)
            gene.add_transcript(transcript)
            pass
        elif line.split('\t')[1] == 'exon':
            id = line.split('\t')[5]
            refstart = int(line.split('\t')[2].split(':')[0])
            refend = int(line.split('\t')[2].split(':')[1])
            start = int(line.split('\t')[3].split(':')[0])
            end = int(line.split('\t')[3].split(':')[1])
            exon = Exon(id, refstart, refend, start, end, transcript)
            transcript.add_exon(exon)
            pass
        line = gtf.readline()
    for chr, genes in chr_to_genes.items():
        tree = IntervalTree()
        for gene in genes:
            for transcript in gene.transcripts:
                start = transcript.start
                end = transcript.end
                tree[start:end] = transcript
        index[chr] = tree
'''
step 1.2
fuction: read_sam_file()
input: thread, ref_index, readfile
output: NULL (global reads)
Use minimap2 to align reads to rebuild transcripts and then index the result of alignment.
'''
def read_sam_file(thread, ref_index, readfile,  output_sam, if_output, type_alignment):
    global reads, overlap_present
    minimap2_view_process = subprocess_popen(shlex.split('minimap2 -ax map-{} -t {} {} {} --secondary=no'.format(type_alignment, str(thread), ref_index, readfile)))
    old_id = ''
    n = 0
    with open(output_sam, 'w') as samfile:
        for line in minimap2_view_process.stdout:
            if if_output == True:
                samfile.write(line)
            if 'SA:' in line and line[0] != '@':
                read_id = line.split('\t')[0]
                read_flag = int(line.split('\t')[1])
                if read_flag != 4 and read_flag != 256 and read_flag != 272:
                    read_chr = line.split('\t')[2]
                    read_start = int(line.split('\t')[3])
                    read_cigar = line.split('\t')[5]
                    read_fasta = line.split('\t')[9]
                    if read_id != old_id:
                        old_id = read_id
                        n += 1
                        if n > 1 and len(read.samlines) > 1:
                            if read_to_transcript(read, overlap_present) == True:
                                reads.append(read)
                        read = Read(read_id, read_start, read_cigar, read_flag, read_chr, read_fasta)
                    else:
                        read.add_samline(read_id, read_start, read_cigar, read_flag, read_chr, read_fasta)
    if len(read.samlines) > 1:
        if read_to_transcript(read, overlap_present) == True:
            reads.append(read)
    minimap2_view_process.stdout.close()
    minimap2_view_process.wait()

def read_sam_file2(output_sam):
    global reads, overlap_present
    old_id = ''
    n = 0
    samfile = open(output_sam)
    line = samfile.readline()
    while line:
        if 'SA:' in line and line[0] != '@':
            read_id = line.split('\t')[0]
            read_flag = int(line.split('\t')[1])
            if read_flag != 4 and read_flag != 256 and read_flag != 272:
                read_chr = line.split('\t')[2]
                read_start = int(line.split('\t')[3])
                read_cigar = line.split('\t')[5]
                read_fasta = line.split('\t')[9]
                if read_id != old_id:
                    old_id = read_id
                    n += 1
                    if n > 1 and len(read.samlines) > 1: #filter3
                        if read_to_transcript(read, overlap_present) == True:
                            reads.append(read)
                    read = Read(read_id, read_start, read_cigar, read_flag, read_chr, read_fasta)
                else:
                    read.add_samline(read_id, read_start, read_cigar, read_flag, read_chr, read_fasta)
        line = samfile.readline()
    if len(read.samlines) > 1:
        if read_to_transcript(read, overlap_present) == True:
            reads.append(read)
    pass

def read_to_transcript(read, precent):
    global index
    for samline in read.samlines:
        Sort = sorted(index[samline.chr][samline.start:samline.end])
        for m in Sort:
            transcriptlength = m[1] - m[0] + 1
            cigarlength = samline.end - samline.start + 1
            if transcriptlength < cigarlength:
                limit = int(transcriptlength * precent)
            else:
                limit = int(cigarlength * precent)
            if samline.end - m[0] >= limit and m[1] - samline.start >= limit:
                samline.add_transcript(m[2], m[0], m[1])
    
    if read.transcripts_and_genes() == True:
        read.num_samline()
        if len(read.genes) > 1:
            return True
    return False

'''
step 1.3
fuction: gene_groups()
input: read_in_same_chr
output: (gene_groups_with_one_read, gene_groups_with_multiple_read)
Find gene pairs with the number of supporting reads larger than one. 
'''               
def gene_groups():
    global reads
    time1 = time.time()
    m = 0
    print('Sparse matrix multiplication-based candidate gene pairing')
    gene_groups = []
    gene_to_vector = {} #{key = gene : value = [read, ...]}
    gene_to_reads = {}
    num_of_reads = len(reads)
    vector_list = []
    for i in range(num_of_reads):
        read = reads[i]
        for gene in list(read.genes.keys()):
            if gene in gene_to_vector:
                gene_to_vector[gene][i] = 1
                gene_to_reads[gene].append(read)
            else:
                gene_to_vector[gene] = np.zeros(num_of_reads, dtype= np.uint8)
                gene_to_vector[gene][i] = 1
                gene_to_reads[gene] = [read,]
    num_of_genes = len(gene_to_vector.keys())
    for gene, vector in gene_to_vector.items():
        vector_list.append(vector)
    gene_Matrix = np.vstack(vector_list)
    csr_gene_Matrix = csr_matrix(gene_Matrix)
    csr_gene_Matrix_T = csr_matrix(gene_Matrix.T)
    groups = csr_gene_Matrix.dot(csr_gene_Matrix_T)
    count = groups.indptr
    index = groups.indices
    data = groups.data
    for i in range(num_of_genes - 1):
        index0 = index[count[i] : count[i + 1]]
        data0 = data[count[i] : count[i + 1]]
        if index0.size > 0:
            for j in range(count[i + 1] - count[i]):
                if data0[j] > 1 and index0[j] > i:
                    gene1 = list(gene_to_vector.keys())[i]
                    gene2 = list(gene_to_vector.keys())[index0[j]]
                    read_list = []
                    add_vector = gene_to_vector[gene1] + gene_to_vector[gene2]
                    index_read = np.where(add_vector == 2)
                    for k in index_read[0]:
                        read_list.append(reads[k])
                    genegroup = Genes_group(gene1= gene1, gene2= gene2, reads= read_list)
                    gene_groups.append(genegroup)
                    m += 1
    time2 = time.time()
    print('used ', time2 - time1, 's')
    return gene_groups
'''
step 1.output_result
fuction: output_result()
input: genegroups, output_file_dir
output: NULL
Output the gene groups to dir.
'''
def output_result(genegroups, output_file_dir):
    time1 = time.time()
    print('-' * 50)
    print('Outputing result')
    if not os.path.exists(output_file_dir):
        os.makedirs(output_file_dir)
    for group in genegroups:
        file_name = group.gene1.name + '::' + group.gene2.name + '.txt'
        with open(output_file_dir + file_name, 'w') as outputfile:
            outputfile.write('@' + group.gene1.id + '::' + group.gene2.id + '\treads num: ' +str(len(group.reads)) + '\n')
            for read in group.reads:
                outputfile.write('>' + read.id + '\n')
                line0 = 'chr\tgene id\ttranscript id\ttstart:end\tstart:end\tposition in read\n'
                length1 = 0
                length2 = 0
                for samline in read.samlines:
                    if samline.transcript.gene.id == group.gene1.id:
                        if samline.length > length1:
                            length1 = samline.length
                            line1 = samline.chr + '\t' + samline.transcript.gene.id + '\t' + samline.transcript.id + '\t' + \
                                str(samline.tstart) + ':' + str(samline.tend) + '\t' + \
                                str(samline.start) + ':' + str(samline.end) + '\t' + str(samline.num) + '\t' + samline.sequence + '\n'
                    elif samline.transcript.gene.id == group.gene2.id:
                        if samline.length > length2:
                            length2 = samline.length
                            line2 = samline.chr + '\t' + samline.transcript.gene.id + '\t' + samline.transcript.id + '\t' + \
                                str(samline.tstart) + ':' + str(samline.tend) + '\t' + \
                                str(samline.start) + ':' + str(samline.end) + '\t' + str(samline.num) + '\t' + samline.sequence + '\n'
                outputfile.write(line0)
                outputfile.write(line1)
                outputfile.write(line2)
    time2 = time.time()
    print('used ', time2 - time1, 's')
'''
step 1.4
fuction: clustering
input: genegroups
output: NULL
Cluster the reads' breakpoints in one gene group.
'''
def clustering(genegroups, tree_dir, out_put, result_out, if_output):
    time1 = time.time()
    print('-' * 50)
    print('Clustering the breakpoint and POA')
    global clustering_num
    nodes0 = []
    n_list = []
    result_file = ''
    POA_file = ''
    with open(out_put, 'w') as output:
        for group in genegroups:
            group.clustering()
            if group.flag == True:
                for flag, tree in group.trees.items():
                    if if_output == True:
                        file_name = group.gene1.name + '::' + group.gene2.name + '_' + str(tree.flag) + '.txt'
                        if not os.path.exists(tree_dir):
                            os.makedirs(tree_dir)
                        tree.print_tree(tree_dir + file_name)
                    nodelist = tree.classify(clustering_num)
                    n = 0
                    for node in nodelist:
                        node.add_group(group, flag)
                        nodes0.append(node)
                        read_list = []
                        read_list_1 = []
                        read_list_2 = []
                        
                        leaf = node.leaves[0]
                        if leaf.anti == True:
                            min_f = len(leaf.bsam.sequence)
                            min_b = len(leaf.fsam.sequence)
                        else:
                            min_b = len(leaf.bsam.sequence)
                            min_f = len(leaf.fsam.sequence)
                        
                        for leaf in node.leaves:
                            if leaf.anti == True:
                                if len(leaf.bsam.sequence) < min_f:
                                    min_f = len(leaf.bsam.sequence)
                                if len(leaf.fsam.sequence) < min_b:
                                    min_b = len(leaf.fsam.sequence)
                                read_list_1.append(Antisense(leaf.bsam.sequence))
                                read_list_2.append(Antisense(leaf.fsam.sequence))
                            else:
                                if len(leaf.fsam.sequence) < min_f:
                                    min_f = len(leaf.fsam.sequence)
                                if len(leaf.bsam.sequence) < min_b:
                                    min_b = len(leaf.bsam.sequence)
                                read_list_1.append(leaf.fsam.sequence)
                                read_list_2.append(leaf.bsam.sequence)
                        min_f_seq = max(min_abpoa, min_f)
                        min_b_seq = max(min_abpoa, min_b)
                        for i in range(len(read_list_1)):
                            seq1 = read_list_1[i]
                            seq2 = read_list_2[i]
                            if min_f_seq < len(seq1):
                                seq1 = seq1[-min_f_seq:]
                            if min_b_seq < len(seq2):
                                seq2 = seq2[:min_b_seq]
                            read_list.append(seq1 + seq2)
                        a = pa.msa_aligner()
                        res=a.msa(read_list, out_cons=True, out_msa=True)
                        for seq in res.cons_seq:
                            line1 = seq + '\n'
                        n += 1
                        line0 = '>' + group.gene1.id + '_' + group.gene2.id + '_' + str(n) + '_' + str(len(node.leaves)) + '\n'
                        output.write(line0)
                        output.write(line1)
                        POA_file += line0 + line1
                        n_list.append(n)
    time2 = time.time()
    print('find ', str(len(nodes0)), 'fusion genes groups')
    print('used ', time2 - time1, 's')
    time1 = time.time()
    print('-' * 50)
    print('Output result')
    groups = {}
    with open(result_out, 'w') as file:
        if if_output == True:
            file.write('gene1 name\tgene2 name\tgene1 id\tgene2 id\tchr1\tchr2\tgene1 pos\tgene2 pos\tbreakpoint1\tbreakpoint2\tnum\tcode\treads\n')
        i = 0
        for node in nodes0:
            code = str(n_list[i])
            i += 1
            gene1 = node.group.gene1.id
            gene2 = node.group.gene2.id
            n1 = node.group.gene1.name
            n2 = node.group.gene2.name
            pos1 = str(node.group.gene1.rstart) + ':' + str(node.group.gene1.rend)
            pos2 = str(node.group.gene2.rstart) + ':' + str(node.group.gene2.rend)
            chr1 = node.group.gene1.chr + ':' + node.group.gene1.flag
            chr2 = node.group.gene2.chr + ':' + node.group.gene2.flag
            read_num = len(node.leaves)
            break1 = int(node.address[0] + node.group.trees[node.gf].mean1)
            break2 = int(node.address[1] + node.group.trees[node.gf].mean2)

            key = '>' + group.gene1.id + '_' + group.gene2.id + '_' + str(n) + '_' + str(len(node.leaves))
            reads_ls = []
            reads_name = ''
            for leaf in node.leaves:
                reads_ls.append(leaf.read)
                reads_name += leaf.read.id+ ';'

            line = n1 + '\t' + n2 + '\t' + \
                gene1 + '\t' + gene2 + '\t' + \
                chr1 + '\t' + chr2 + '\t' + \
                pos1 + '\t' + pos2 + '\t' + \
                str(break1) + '\t' + str(break2) + '\t' + \
                str(read_num) + '\t' + code + '\t' + reads_name +'\n'
            if if_output == True:
                file.write(line)
            result_file += line
            if key not in groups:
                groups[key] = reads_ls
    time2 = time.time()
    print('used ', time2 - time1, 's')
    print('-' * 50)
    return result_file, POA_file, groups
'''
fuction: usage_and_exit
'''
def usage_and_exit():
    sys.stderr.write('Usage: genefusion [indexdir]')
    sys.stderr.write('\n')
    exit(0)
'''
fuction input_files
'''
def input_files():
    if len(sys.argv) < 3:
        usage_and_exit()
    else:
        annotation = sys.argv[1] + 'rebuild_annotation.txt'
        ref_index = sys.argv[1] + 'trans_index.mmi'
        readfile = sys.argv[2]
        if len(sys.argv) > 3:
            thread = int(sys.argv[3])
        else:
            thread = 2
    return (annotation, ref_index, readfile, thread)

'''
fuction: control_thread
input: annotation, ref_index, readfile, thread, if_output, link_tuple
output: NULL
Make sure the thread number is posible.
'''
def control_thread(annotation, ref_index, readfile, thread, if_output, link_tuple, cnums, if_countine):
    global index
    global reads
    global min_length, exon_boundary, overlap_present, clustering_num, sup_read_num
    output_sam, output_file_dir, tree_dir, out_put, result_out = link_tuple
    min_length, exon_boundary, overlap_present, clustering_num , type_alignment, sup_read_num= cnums
    cpu = os.cpu_count()
    if cpu < thread:
        thread = cpu
    time1 = time.time()
    print('-' * 50)
    print('Reading index file')
    read_annotation(annotation)
    time2 = time.time()
    print('used ', time2 - time1, 's')
    print('-' * 50)
    time1 = time.time()
    print('Transcriptome-based alignment')
    print('-' * 50)
    if if_countine == False:
        read_sam_file(thread, ref_index, readfile, output_sam, if_output, type_alignment)
    else:
        read_sam_file2(output_sam)
    time2 = time.time()
    print('found ', len(reads), ' reads with genefusion')
    print('used ', time2 - time1, 's')
    print('-' * 50)
    genegroups = gene_groups()
    if if_output == True:
        output_result(genegroups, output_file_dir)
    result_file, POA_file, groups = clustering(genegroups, tree_dir, out_put, result_out, if_output)
    return result_file, POA_file, groups
'''
fuction: main
'''
def main(readfile, index_dir, middlefile, thread, if_output, cnums, if_countine):
    annotation = index_dir + 'rebuild_annotation.txt'
    ref_index = index_dir + 'trans_index.mmi'
    output_sam = middlefile + 'step1.DNA_like_alignment.sam'
    output_file_dir = middlefile + 'step1.gene_groups/'
    tree_dir = middlefile + 'step1.clustering_tree/'
    out_put = middlefile + 'step1.POA_result.fasta'
    result_out = middlefile + 'step1.genefusion.csv'
    link_tuple = (output_sam, output_file_dir, tree_dir, out_put, result_out)
    time1 = time.time()
    result_file, POA_file, groups = control_thread(annotation, ref_index, readfile, thread, if_output, link_tuple, cnums, if_countine)
    time2 = time.time()
    timeuse = time2 - time1
    print('Transcriptome-based detection totally used ' + str(timeuse) + 's')
    print('-' * 50)
    return result_file, POA_file, groups
'''
main fuction
'''
if __name__ == "__main__":
    annotation, ref_index, readfile, thread = input_files()
    link_tuple = (output_sam, output_file_dir, tree_dir, out_put, result_out)
    print('-' * 50)
    print('Start detecting!!')
    time1 = time.time()
    cnums = (min_length, exon_boundary, overlap_present, clustering_num, type_alignment)
    control_thread(annotation, ref_index, readfile, thread, True, link_tuple, cnums)
    time2 = time.time()
    timeuse = time2 - time1
    print('Totally using ' + str(timeuse) + 's')
    print('-' * 50)