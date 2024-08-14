import time
import os, sys
import argparse
import step0_preprocessing
import step1_alignment_flitering_clustering
import step3_Recheck_and_Classify

'''
fuction: sub_cmd_index
'''
def sub_cmd_index(args):
    print('#########################')
    print('Start indexing')
    time1 = time.time()
    step0_preprocessing.main(vars(args)['<annotationfile>'], vars(args)['<referencefile>'], vars(args)['<indexdir>'])
    os.system('minimap2 -d ' + vars(args)['<indexdir>'] + 'ref_index.mmi ' + vars(args)['<referencefile>'])
    os.system('minimap2 -d ' + vars(args)['<indexdir>'] + 'trans_index.mmi ' + vars(args)['<indexdir>'] + 'rebuild_transcript.fa')
    time2 = time.time()
    timeuse = time2 - time1
    print('Totally using ' + str(timeuse) + 's')
    print('#########################')
'''
fuction: sub_cmd_index
'''
def sub_cmd_detect(args):
    print('###########################################')
    print('Start detecting!!')
    time1 = time.time()
    if not os.path.exists(args.middlefile):
        os.makedirs(args.middlefile)
    #step1
    print('###########################################')
    print('step1: DNA-like alignment, Filtering SA signal, Clustering the breakpoint and POA')
    step1_alignment_flitering_clustering.main(vars(args)['<readfile>'], vars(args)['<indexdir>'], args.middlefile, args.thread, args.print_middle_output)
    #step2
    print('###########################################')
    print('step2: RNA-seq alignment')
    out_put = args.middlefile + 'step1.POA_result.fasta'
    sam_file3 = args.middlefile + 'step2.RNA_alignment.sam'
    sam_file4 = args.middlefile + 'step2.RNA_alignment.only_primary.sam'
    os.system('minimap2 -ax splice -uf -k14 ' + (vars(args)['<indexdir>'] + 'ref_index.mmi ') + out_put + ' -o ' + sam_file3)
    #os.system('samtools view -h -F 4 '+ sam_file3 +' -o ' + sam_file4)
    #step3
    print('###########################################')
    print('step3: Checking the result')
    #print(args.output)
    step3_Recheck_and_Classify.main(args.output + '.csv', args.middlefile)
    if args.print_middle_output == False:
        os.system('rm -rf ' + args.middlefile)
    time2 = time.time()
    timeuse = time2 - time1
    print('Totally used ' + str(timeuse) + 's')
    print('#########################')
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
    detect_parser.add_argument('-m', '--middlefile', type= str, default= './middlefile/', help= 'Middle file directory (default= "./middlefile/")', metavar= 'dir')
    detect_parser.add_argument('-o', '--output', type= str, default= './result', help= 'Output result name (example: if input "-o result", then GFHunter will output result.csv) (default= "./result")', metavar= 'str')
    detect_parser.add_argument('-M', '--print_middle_output', action= 'store_true', help= 'if you want to save middle file, you can add the option')
    detect_parser.add_argument('-t', '--thread', type= int, default= 4, help= 'Threads GFHunter used', metavar= 'int')
    detect_parser.set_defaults(func= sub_cmd_detect)
    
    args = parser.parse_args()
    #print(vars(args))
    if not hasattr(args, 'func'):
        args = parser.parse_args(['-h'])
    args.func(args)
if __name__ == "__main__":
    Parser_set()