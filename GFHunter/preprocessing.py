import time
import os, sys

rebuild_transcript_file = 'rebuild_transcript.fa'
rebuild_index_file = 'rebuild_annotation.txt'

'''
class: Gene
'''
class Gene():
    def __init__(self, id, name, start, end, chr, flag, gene_type):
        self.id = id
        self.name = name
        self.start = start
        self.end = end
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
    def __init__(self, id, name, start, end, chr, flag, trans_type):
        self.id = id
        self.name = name
        self.start = start
        self.end = end
        self.chr = chr
        self.flag = flag
        self.trans_type = trans_type
        self.exons = []
        pass
    def add_exon(self, exon):
        self.exons.append(exon)
        pass
    def add_pos(self, start, end):
        self.restart = start
        self.reend = end
'''
class: Exon
'''
class Exon():
    def __init__(self, id, start, end):
        self.id = id
        self.start = start
        self.end = end
        pass
    def add_sequence(self, sequence):
        self.sequence = sequence
        pass
    def add_pos(self, start, end):
        self.restart = start
        self.reend = end
'''
fuction: Read_Annotation
input: gtffile
output: chr_to_genes #{chr : [gene1, gene2, ...]}
'''
def Read_Annotation(gtffile):
    print('-' * 50)
    print('Reading annotation file')
    time1 = time.time() 

    m = 0
    n = 0
    chr_to_genes = {}
    gtf = open(gtffile)
    line = gtf.readline()
    while line:
        if '##' not in line:
            if line.split('\t')[2] == 'gene':
                m += 1
                id = line.split('\t')[8].split(';')[0].split('"')[1]
                name = line.split('\t')[8].split(';')[2].split('"')[1]
                start = int(line.split('\t')[3])
                end = int(line.split('\t')[4])
                chr = line.split('\t')[0]
                flag = line.split('\t')[6]
                gene_type = line.split('\t')[8].split(';')[1].split('"')[1]
                gene = Gene(id, name, start, end, chr, flag, gene_type)
                if chr in chr_to_genes:
                    chr_to_genes[chr].append(gene)
                else:
                    chr_to_genes[chr] = [gene,]
                pass
            elif line.split('\t')[2] ==  'transcript':
                n += 1
                id = line.split('\t')[8].split(';')[1].split('"')[1]
                name = line.split('\t')[8].split(';')[5].split('"')[1]
                start = int(line.split('\t')[3])
                end = int(line.split('\t')[4])
                chr = line.split('\t')[0]
                flag = line.split('\t')[6]
                trans_type = line.split('\t')[8].split(';')[4].split('"')[1]
                transcript = Transcript(id, name,start, end, chr, flag, trans_type)
                gene.add_transcript(transcript)
                pass
            elif line.split('\t')[2] == 'exon':
                id = line.split('\t')[8].split(';')[6].split(' ')[2]
                start = int(line.split('\t')[3])
                end = int(line.split('\t')[4])
                exon = Exon(id, start, end)
                transcript.add_exon(exon)
                pass
        line = gtf.readline()
    
    time2 = time.time()
    timeuse = time2 - time1
    print('Using ' + str(timeuse) + 's')
    print('Read ', m, ' genes')
    print('Read ', n, ' transcriprts')
    print('-' * 50)
    return chr_to_genes
'''
fuction: Filtering
input: chr_to_genes #{chr : [gene1, gene2, ...]}
output: NULL
'''
def Filtering(chr_to_genes):
    print('-' * 50)
    print('Filtering the genes')
    time1 = time.time() 
    new_chr_to_genes = {}
    m = 0
    n = 0
    for chr, genes in chr_to_genes.items():
        if chr != 'chrM' and chr != 'chrY':
            new_chr_to_genes[chr] = []
            for gene in genes:
                if gene.type == 'protein_coding':
                    new_chr_to_genes[chr].append(gene)
                    m += 1
                    for transcript in gene.transcripts:
                        n += 1
                        if gene.flag == '-':
                            transcript.exons = transcript.exons[::-1]
                        '''for exon in transcript.exons:
                            pass'''
    
    time2 = time.time()
    timeuse = time2 - time1
    print('Using ' + str(timeuse) + 's')
    print('Reserved ', m, ' genes')
    print('Reserved ', n, ' transcriprts')
    print('-' * 50)
    return new_chr_to_genes
'''
fuction: Read_Reference
input: reference, chr_to_genes
output: NULL
'''
def Read_Reference(reference, chr_to_genes):
    print('-' * 50)
    print('Read reference file')
    time1 = time.time() 

    chrs = chr_to_genes.keys()
    ref = open(reference)
    line = ref.readline()
    while line:
        if '>' in line:
            chr = line.split('>')[1].split(' ')[0]
            if chr in chrs:
                sequence = ''
                line = ref.readline()
                while '>' not in line:
                    sequence += line.split('\n')[0]
                    line = ref.readline()
                for gene in chr_to_genes[chr]:
                    for transcript in gene.transcripts:
                        for exon in transcript.exons:
                            exon.add_sequence(sequence[exon.start - 1: exon.end])
            else:
                line = ref.readline()
        else:
            line = ref.readline()
    
    time2 = time.time()
    timeuse = time2 - time1
    print('Using ' + str(timeuse) + 's')
    print('-' * 50)
'''
fuction: Build_transcript_file
input: chr_to_genes, rebuild_transcript_file, rebuild_index_file, middlefile
output: NULL
'''
def Build_transcript_file(chr_to_genes, rebuild_transcript_file, rebuild_index_file, middlefile):
    print('-' * 50)
    print('Building transcript file and index')
    time1 = time.time() 

    if not os.path.exists(middlefile):
        os.makedirs(middlefile)
    with open(middlefile + rebuild_transcript_file, 'w') as rtf:
        for chr, genes in chr_to_genes.items():
            rtf.write('>' + chr + '\n')
            refline = ''
            length = 0
            for gene in genes:
                for transcript in gene.transcripts:
                    tstart = length + 1
                    for exon in transcript.exons:
                        start = length + 1
                        refline += exon.sequence
                        length += exon.end - exon.start + 1
                        end = length
                        exon.add_pos(start, end)
                    tend = length
                    transcript.add_pos(tstart, tend)
            one_line = ''
            for i in refline:
                one_line += i
                if len(one_line) == 100:
                    rtf.write(one_line + '\n')
                    one_line = ''
            if one_line != '':
                rtf.write(one_line + '\n')
    with open(middlefile + rebuild_index_file, 'w') as rif:
        line0 = 'chr\tclass\trefpos\trebuildpos\tflag\tid\tname\ttype\n'
        rif.write(line0)
        for chr, genes in chr_to_genes.items():
            length = 0
            for gene in genes:
                line = gene.chr + '\tgene\t' + \
                    str(gene.start) + ':' + str(gene.end) + \
                    '\t-\t' + \
                    gene.flag + '\t' + gene.id + '\t' + gene.name + '\t' + \
                    gene.type + '\n'
                rif.write(line)
                for transcript in gene.transcripts:
                    line = transcript.chr + '\ttranscript\t' + \
                        str(transcript.start) + ':' + str(transcript.end) + '\t' + \
                        str(transcript.restart) + ':' + str(transcript.reend) + '\t' + \
                        transcript.flag + '\t' + transcript.id + '\t' + transcript.name + '\t' + \
                        transcript.trans_type+'\n'
                    rif.write(line)
                    for exon in transcript.exons:
                        line = transcript.chr + '\texon\t' + \
                            str(exon.start) + ':' + str(exon.end) + '\t' + \
                            str(exon.restart) + ':' + str(exon.reend) + '\t' + \
                            transcript.flag + '\t' + exon.id +'\n'
                        rif.write(line)
    
    time2 = time.time()
    timeuse = time2 - time1
    print('Using ' + str(timeuse) + 's')
    print('-' * 50)
'''
fuction: usage_and_exit
'''
def usage_and_exit():
    sys.stderr.write('Usage: index [annotationfile] [referencefile] [indexdir]')
    sys.stderr.write('\n')
    exit(0)
'''
fuction: input_files
'''
def input_files():
    if len(sys.argv) < 4:
        usage_and_exit()
    else:
        gtffile = sys.argv[1]
        reference = sys.argv[2]
        middlefile = sys.argv[3]
    return (gtffile, reference, middlefile)
'''
fuction: main
'''
def main(gtffile, reference, middlefile):
    chr_to_genes = Read_Annotation(gtffile)
    chr_to_genes = Filtering(chr_to_genes)
    Read_Reference(reference, chr_to_genes)
    Build_transcript_file(chr_to_genes, rebuild_transcript_file, rebuild_index_file, middlefile)
'''
main fuction
'''
if __name__ == "__main__":
    gtffile, reference, middlefile = input_files()
    chr_to_genes = Read_Annotation(gtffile)
    chr_to_genes = Filtering(chr_to_genes)
    Read_Reference(reference, chr_to_genes)
    Build_transcript_file(chr_to_genes, rebuild_transcript_file, rebuild_index_file, middlefile)