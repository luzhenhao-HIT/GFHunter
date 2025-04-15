# GFHunter enables accurate and efficient gene fusion detection in long-read cancer transcriptomes
## Overview
- [Introduction](#introduction)
- [Dependence](#dependence)
- [Installation](#installation)
- [Usage](#usage)
- [Citation](#citation)
- [Contact](#contact)

## Introduction
- **GFHunter** is a long read sequence transcriptome alignment-based fusion genes detection tool.
- **Indexing**: User can download the prepared index folder based on **GENCODE release 44** here or create a new index by GFHunter **index** function.
- **Detection**: GFHunter detects fusions in **long read RNA-seq data**, user just need **index folder** and **long read data** under detection.
- Benchmark evaluations position GFHunter as redefining the analytical frontiers of gene fusion detection:
  1. **Lowest false-positive fusion rates** in non-tumor data.
  2. Achieves a **substantially higher F1-score** in real cancer cell lines evaluation.
  3. Operates **many times faster** and requires **fewer and stable memory** than other tools.
## Dependence
```
1. python >=3
2. minimap2 >=2.22
3. intervaltree =3.1.0
4. numpy >=2.0.1
5. scipy >=1.14.0
6. pyabpoa =1.5.2
7. gcc >=6.4.0
8. pysam >= 0.23.0
```
## Installation
```
#install minimap2

conda create -n GFHunter python=3.8
conda activate GFHunter
conda install bioconda::minimap2

#install GFHunter

pip install GFHunter
```
## Usage
GFHunter offers 2 steps to detect fusions: **index** and **detect**
```
usage: GFHunter [-h] {index,detect,sc} ...

positional arguments:
  {index,detect,sc}
    index            Create index of the GFHunter
    detect           Dectect gene fusions
    sc               Dectect gene fusions on single cell RNA-seq data

optional arguments:
  -h, --help         show this help message and exit
```
### Index
```
GFHunter index [-h] <annotationfile> <referencefile> <indexdir>

positional arguments:
  <annotationfile>  Gene annotation file (gtf)
  <referencefile>   Reference file (fasta/fa)
  <indexdir>        Output index directory

options:
  -h, --help        show this help message and exit
```
Note: the annotation only support GENCODE currently.
### Detection
```
usage: GFHunter detect [-h] [-o str] [-m dir] [-M] [-t int] [-T type] [-n int] [-e int] [-p float] [-c int] [-l int] [-b int] [-C] <readfile> <indexdir>

positional arguments:
  <readfile>            Read file (fasta/fastq)
  <indexdir>            Index directory

optional arguments:
  -h, --help            show this help message and exit
  -o str, --output str  the name of fusion result (default = "./result")
  -m dir, --middlefile dir
                        temporary folder of middle files (default = "./middlefile/")
  -M, --print_middle_output
                        retain the middle files in detection
  -t int, --threads int
                        threads GFHunter used (default = 4)
  -T type, --trans_based_align_type type
                        setting the transcriptome-based alignmnet type of minimap2: pb/hifi/ont/iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping (default = ont)
  -n int, --min_read_length int
                        minimum length of reads (bp) considered (default = 50)
  -e int, --max_exon_boundary int
                        maximum length between breakpoint and exon boundary (default = 30)
  -p float, --overlap_precent float
                        precent of overlap between reads and transcripts (default = 0.5)
  -c int, --min_clustering_length int
                        minimum distance between two cluster (default = 200)
  -l int, --least_support_reads int
                        least reads number to support gene fusions (default = 2)
  -b int, --max_breakpoint_distance int
                        maximum distance between crossing refinement breakpoints and transcriptome based detection breakpoints (default = 200)
  -C, --countine        countine detection)
```
| Parameter | Description | Default |
|-----|---------------|-----|
|--output|the name of fusion result|./result|
|--middlefile|temporary folder of middle files|./middlefile/|
|--print_middle_output|retain the middle files in detection|NULL|
|--thread|threads GFHunter used|4|
|--trans_based_align_type|setting the transcriptome-based alignmnet type of minimap2: pb/hifi/ont/iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping|ont|
|--min_read_length|minimum length of reads (bp) considered|50|
|--max_exon_boundary|maximum length between breakpoint and exon boundary|30|
|--overlap_precent|precent of overlap between reads and transcripts|0.5|
|--min_clustering_length|minimum distance between two cluster|200|
|--least_support_reads|least reads number to support gene fusions|2|
### Single-cell detection
If user want to detect fusions on single cell RNA-seq data, they should use wf-single-cell to preprocess the raw fastq. After the preprocess, use the tagged bam as input of GFHunter with the fuction "sc".
```
usage: GFHunter sc [-h] [-o str] [-m dir] [-M] [-t int] [-T type] [-n int] [-e int] [-p float] [-c int] [-l int] [-L int] [-b int] [-C] <bamfile> <indexdir>

positional arguments:
  <bamfile>             Read file handled by wf-single-cell (bam)
  <indexdir>            Index directory

optional arguments:
  -h, --help            show this help message and exit
  -o str, --output str  the name of fusion result (default = "./result")
  -m dir, --middlefile dir
                        temporary folder of middle files (default = "./middlefile/")
  -M, --print_middle_output
                        retain the middle files in scion
  -t int, --threads int
                        threads GFHunter used (default = 4)
  -T type, --trans_based_align_type type
                        setting the transcriptome-based alignmnet type of minimap2: pb/hifi/ont/iclr - CLR/HiFi/Nanopore/ICLR vs reference mapping (default = ont)
  -n int, --min_read_length int
                        minimum length of reads (bp) considered (default = 50)
  -e int, --max_exon_boundary int
                        maximum length between breakpoint and exon boundary (default = 30)
  -p float, --overlap_precent float
                        precent of overlap between reads and transcripts (default = 0.5)
  -c int, --min_clustering_length int
                        minimum distance between two cluster (default = 200)
  -l int, --least_support_reads int
                        least reads number to support gene fusions (default = 2)
  -L int, --least_support_cells int
                        least cells number to support gene fusions (default = 5)
  -b int, --max_breakpoint_distance int
                        maximum distance between crossing refinement breakpoints and transcriptome based scion breakpoints (default = 200)
  -C, --countine        countine detection
```
### Note
- GFHunter only support **GENCODE annotation** at present.
- User can dowload a pre-prepared index from https://drive.google.com/file/d/1UuWu0bVQm7lUh7AxVyQDC5zFQe37iOHv/view?usp=drive_link.
## Citation
## Contact
