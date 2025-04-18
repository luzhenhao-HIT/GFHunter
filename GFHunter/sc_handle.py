import pysam
import os
import time

class Cell():
    def __init__(self, CB):
        self.CB = CB
        self.reads = []
        self.fusions = []
        pass
    def add_reads(self, read):
        self.reads.append(read)
    def add_fusions(self, fusion):
        if fusion not in self.fusions:
            self.fusions.append(fusion)

def read_bam(outfile, bam, output):
    time1 = time.time()
    f_name, e = os.path.splitext(os.path.basename(bam))
    indexfile = outfile + f_name + '.CBindex'
    fasta = output + '.fasta'
    bamfile = pysam.AlignmentFile(bam, 'rb')
    anti_index = {}
    if not os.path.exists(indexfile):
        with open(fasta, 'w') as fasta_out:
            for read in bamfile:
                CB = read.get_tag('CB')
                UMI = read.get_tag('UB')
                if CB in anti_index:
                    anti_index[CB].append((read.query_name, UMI))
                else:
                    anti_index[CB] = [(read.query_name, UMI),]
                if read.flag in [0, 16]:
                    fasta_out.write(f">{read.query_name}\n{read.query_sequence}\n")
                

        with open(indexfile, 'w') as f:
            f.write('CB\tread\n')
            for CB, reads in anti_index.items():
                for read in reads:
                    name, UMI = read
                    line = CB + '\t' + name + '\t' + UMI + '\n'
                    f.write(line)
    else:
        print('Index exists, jump to the next step')
    time2 = time.time()
    timeuse = time2 - time1
    print('Indexing time used: ', timeuse)

def create_key(a, b):
    ls = [a, b]
    ls.sort()
    key = ls[0] + ':' + ls[1]
    return key

def read_index(file):
    index = {}
    cells = {}
    f = open(file)
    line = f.readline()
    while line:
        CB = line.split('\t')[0]
        read = line.split('\t')[1].split('\n')[0]
        if CB in cells:
            cells[CB].add_reads(read)
        else:
            cell = Cell(CB)
            cells[CB] = cell
        index[read] = cells[CB]
        line = f.readline()
    return index

def read_fusions(file, index, output, lc):
    fusion_cells = []
    fusion_line = {}
    result = open(file)
    line = result.readline()
    line = result.readline()
    while line:
        g1 = line.split('\t')[0]
        g2 = line.split('\t')[1]
        id1 = line.split('\t')[2]
        id2 = line.split('\t')[3]
        bps = line.split('\t')[4]
        gp1 = line.split('\t')[5]
        gp2 = line.split('\t')[6]
        rnum = line.split('\t')[7]
        flag = line.split('\t')[8]
        reads_l = line.split('\t')[9]
        key = create_key(g1, g2)
        reads = line.split('\t')[9].split('\n')[0].split(';')[:-1]
        cells = []
        for read in reads:
            cell = index[read]
            cells.append(cell.CB)
            if flag == 'flag = 5: Reliable Fusion':
                cell.add_fusions(key)
                if cell not in fusion_cells:
                    fusion_cells.append(cell)
        cs = list(set(cells))
        cell_l = ''
        for c in cs:
            cell_l += c + ';'
        new_line = g1 + '\t' + g2 + '\t' + id1 + '\t' + id2 + '\t' + bps + '\t' + gp1 + '\t' + gp2 + '\t' + str(len(cs)) + '\t' + rnum + '\t' + flag + '\t' + cell_l + '\t' + reads_l
        if key not in fusion_line:
            fusion_line[key] = [new_line,]
        else:
            fusion_line[key].append(new_line)
        line = result.readline()

    with open(output + '.sc.tsv', 'w') as outfile:
        outfile.write('gene1 name\tgene2 name\tgene1 id\tgene2 id\tbreak points\tgene1 position\tgene2 position\t#cells\t#reads\trank class\tcells\treads\n')
        for fusion, lines in fusion_line.items():
            sum_c_num = 0
            for line in lines:
                c_num = int(line.split('\t')[7])
                sum_c_num += c_num
            if sum_c_num >= lc:
                for line in lines:
                    outfile.write(line)
    return fusion_cells

def output_result(fusion_cells, file, lc):
    index = {}
    with open(file + '.r1.tsv', 'w') as rs:
        rs.write('cell\tfusions\n')
        for cell in fusion_cells:
            l = ''
            for f in cell.fusions:
                l += f + ';'
                if f in index:
                    index[f].append(cell.CB)
                else:
                    index[f] = [cell.CB,]
            line = cell.CB + '\t' + l + '\n'
            rs.write(line)
    sorted_index = dict(sorted(index.items(), key=lambda item: len(item[1]), reverse=True))
    with open(file + '.r2.tsv', 'w') as rs:
        rs.write('fusion\tcell\n')
        for f, cells in sorted_index.items():
            if len(cells) < lc:
                break
            for c in cells:
                line = f + '\t' + c + '\n'
                rs.write(line)
    with open(file + '.sum.tsv', 'w') as rs:
        rs.write('fusion\t#cell\tcells\n')
        for f, cells in sorted_index.items():
            if len(cells) < lc:
                break
            l = ''
            for c in cells:
                l += c + ';'
            line = f + '\t' + str(len(cells)) + '\t' + l + '\n'
            rs.write(line)

def handle(indexfile, bam, output, lc):
    f_name, e = os.path.splitext(os.path.basename(bam))
    indexf = indexfile + f_name + '.CBindex'
    result = output + '.tsv'
    index = read_index(indexf)
    fusion_cells = read_fusions(result, index, output, lc)
    output_result(fusion_cells, output, lc)