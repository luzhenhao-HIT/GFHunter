import time
import os
import argparse
from . import preprocessing
from . import transcriptome_based_detection
from . import crossing_reference_based_refinement
from . import sc_handle


'''
fuction: sub_cmd_index
'''
def sub_cmd_index(args):
    time1 = time.time()
    preprocessing.main(vars(args)['<annotationfile>'], vars(args)['<referencefile>'], vars(args)['<indexdir>'])
    os.system('minimap2 -d ' + vars(args)['<indexdir>'] + 'ref_index.mmi ' + vars(args)['<referencefile>'])
    os.system('minimap2 -d ' + vars(args)['<indexdir>'] + 'trans_index.mmi ' + vars(args)['<indexdir>'] + 'rebuild_transcript.fa')
    time2 = time.time()
    timeuse = time2 - time1
    print('Totally using ' + str(timeuse) + 's')
    print('-' * 50)
'''
fuction: sub_cmd_detect
'''
def sub_cmd_detect(args):
    time1 = time.time()
    if not os.path.exists(args.middlefile):
        os.makedirs(args.middlefile)
    #step1
    print('-' * 50)
    print('TRANSCRIPTOME-BASED DETECTION')
    cnums = (args.min_read_length, args.max_exon_boundary, args.overlap_precent, args.min_clustering_length, args.trans_based_align_type, args.least_support_reads)
    result_file, POA_file, supporting_reads = transcriptome_based_detection.main(vars(args)['<readfile>'], vars(args)['<indexdir>'], args.middlefile, args.threads, args.print_middle_output, cnums, args.countine)
    #step2
    print('-' * 50)
    print('CROSS-REFERENCING-BASED FUSION REFINEMENT')
    print('-' * 50)
    POA_file = args.middlefile + 'step1.POA_result.fasta'
    crossing_reference_based_refinement.main(args.output + '.csv', POA_file, result_file, vars(args)['<indexdir>'], args.middlefile, args.threads, args.print_middle_output, args.least_support_reads, args.max_breakpoint_distance)
    os.system('mv ' + POA_file + ' ' + args.output + '.fusion_transcripts.fasta')
    if args.print_middle_output == False:
        os.system('rm -rf ' + args.middlefile)
    time2 = time.time()
    timeuse = time2 - time1
    print('-' * 50)
    print('Totally used ' + str(timeuse) + 's')
    print('-' * 50)
'''
fuction: sub_cmd_sc
'''
def sub_cmd_sc(args):
    time1 = time.time()
    if not os.path.exists(args.middlefile):
        os.makedirs(args.middlefile)
    #step0
    print('-' * 50)
    print('SINGLE CELL DATA INDEXING')
    sc_handle.read_bam(vars(args)['<indexdir>'], vars(args)['<bamfile>'], args.output)
    #step1
    print('-' * 50)
    print('TRANSCRIPTOME-BASED DETECTION')
    cnums = (args.min_read_length, args.max_exon_boundary, args.overlap_precent, args.min_clustering_length, args.trans_based_align_type, args.least_support_reads)
    result_file, POA_file, supporting_reads = transcriptome_based_detection.main(args.output + '.fasta', vars(args)['<indexdir>'], args.middlefile, args.threads, args.print_middle_output, cnums, args.countine)
    #step2
    print('-' * 50)
    print('CROSS-REFERENCING-BASED FUSION REFINEMENT')
    print('-' * 50)
    POA_file = args.middlefile + 'step1.POA_result.fasta'
    crossing_reference_based_refinement.main(args.output + '.tsv', POA_file, result_file, vars(args)['<indexdir>'], args.middlefile, args.threads, args.print_middle_output, args.least_support_reads, args.max_breakpoint_distance)
    os.system('mv ' + POA_file + ' ' + args.output + '.fusion_transcripts.fasta')
    if args.print_middle_output == False:
        os.system('rm -rf ' + args.middlefile)
    #step3
    print('-' * 50)
    print('MAP FUSION TO CELLS')
    print('-' * 50)
    sc_handle.handle(vars(args)['<indexdir>'], vars(args)['<bamfile>'], args.output, args.least_support_cells)
    time2 = time.time()
    timeuse = time2 - time1
    print('-' * 50)
    print('Totally used ' + str(timeuse) + 's')
    print('-' * 50)

'''
fuction: 
'''
def Parser_set():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    index_parser = subparser.add_parser('index', help = 'Create index of the GFHunter')
    index_parser.add_argument('<annotationfile>', type= str, help= 'Gene annotation file (gtf)')
    index_parser.add_argument('<referencefile>', type= str, help= 'Reference file (fasta/fa)')
    index_parser.add_argument('<indexdir>', type= str, default='./index/', help= 'Output index directory')
    index_parser.set_defaults(func= sub_cmd_index)

    detect_parser = subparser.add_parser('detect', help = 'Dectect gene fusions')
    detect_parser.add_argument('<readfile>', type= str, help= 'Read file (fasta/fastq)')
    detect_parser.add_argument('<indexdir>', type= str, default='./index/', help= 'Index directory')
    detect_parser.add_argument('-o', '--output', type= str, default= './result', help= 'the name of fusion result (default = "./result")', metavar= 'str')
    detect_parser.add_argument('-m', '--middlefile', type= str, default= './middlefile/', help= 'temporary folder of middle files (default = "./middlefile/")', metavar= 'dir')
    detect_parser.add_argument('-M', '--print_middle_output', action= 'store_true', help= 'retain the middle files in detection')
    detect_parser.add_argument('-t', '--threads', type= int, default= 4, help= 'threads GFHunter used (default = 4)', metavar= 'int')
    detect_parser.add_argument('-T', '--trans_based_align_type', type= str, default= 'ont', help= 'setting the transcriptome-based alignmnet type of minimap2: pb/hifi/ont/iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping (default = ont)', metavar= 'type')
    detect_parser.add_argument('-n', '--min_read_length', type= int, default= 50, help= 'minimum length of reads (bp) considered (default = 50)', metavar= 'int')
    detect_parser.add_argument('-e', '--max_exon_boundary', type= int, default= 30, help= 'maximum length between breakpoint and exon boundary (default = 30)', metavar= 'int')
    detect_parser.add_argument('-p', '--overlap_precent', type= float, default= 0.5, help= 'precent of overlap between reads and transcripts (default = 0.5)', metavar= 'float')
    detect_parser.add_argument('-c', '--min_clustering_length', type= int, default= 200, help= 'minimum distance between two cluster (default = 200)', metavar= 'int')
    detect_parser.add_argument('-l', '--least_support_reads', type= int, default= 2, help= 'least reads number to support gene fusions (default = 2)', metavar= 'int')
    detect_parser.add_argument('-b', '--max_breakpoint_distance', type= int, default= 200, help= 'maximum distance between crossing refinement breakpoints and transcriptome based detection breakpoints (default = 200)', metavar= 'int')
    detect_parser.add_argument('-C', '--countine', action= 'store_true', help= 'countine detection')

    detect_parser.set_defaults(func= sub_cmd_detect)

    sc_parser = subparser.add_parser('sc', help = 'Dectect gene fusions on single cell RNA-seq data')
    sc_parser.add_argument('<bamfile>', type= str, help= 'Read file handled by wf-single-cell (bam)')
    sc_parser.add_argument('<indexdir>', type= str, default='./index/', help= 'Index directory')
    sc_parser.add_argument('-o', '--output', type= str, default= './result', help= 'the name of fusion result (default = "./result")', metavar= 'str')
    sc_parser.add_argument('-m', '--middlefile', type= str, default= './middlefile/', help= 'temporary folder of middle files (default = "./middlefile/")', metavar= 'dir')
    sc_parser.add_argument('-M', '--print_middle_output', action= 'store_true', help= 'retain the middle files in scion')
    sc_parser.add_argument('-t', '--threads', type= int, default= 4, help= 'threads GFHunter used (default = 4)', metavar= 'int')
    sc_parser.add_argument('-T', '--trans_based_align_type', type= str, default= 'ont', help= 'setting the transcriptome-based alignmnet type of minimap2: pb/hifi/ont/iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping (default = ont)', metavar= 'type')
    sc_parser.add_argument('-n', '--min_read_length', type= int, default= 50, help= 'minimum length of reads (bp) considered (default = 50)', metavar= 'int')
    sc_parser.add_argument('-e', '--max_exon_boundary', type= int, default= 30, help= 'maximum length between breakpoint and exon boundary (default = 30)', metavar= 'int')
    sc_parser.add_argument('-p', '--overlap_precent', type= float, default= 0.5, help= 'precent of overlap between reads and transcripts (default = 0.5)', metavar= 'float')
    sc_parser.add_argument('-c', '--min_clustering_length', type= int, default= 200, help= 'minimum distance between two cluster (default = 200)', metavar= 'int')
    sc_parser.add_argument('-l', '--least_support_reads', type= int, default= 2, help= 'least reads number to support gene fusions (default = 2)', metavar= 'int')
    sc_parser.add_argument('-L', '--least_support_cells', type= int, default= 5, help= 'least cells number to support gene fusions (default = 5)', metavar= 'int')
    sc_parser.add_argument('-b', '--max_breakpoint_distance', type= int, default= 200, help= 'maximum distance between crossing refinement breakpoints and transcriptome based scion breakpoints (default = 200)', metavar= 'int')
    sc_parser.add_argument('-C', '--countine', action= 'store_true', help= 'countine detection')

    sc_parser.set_defaults(func= sub_cmd_sc)
    
    args = parser.parse_args()
    if not hasattr(args, 'func'):
        args = parser.parse_args(['-h'])
    args.func(args)
if __name__ == "__main__":
    Parser_set()