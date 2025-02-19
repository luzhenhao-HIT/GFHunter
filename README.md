# GFHunter: A Computational Framework for Precision Detection of Gene Fusions in Long-Read Cancer Transcriptomes
## Overview

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
1. python>=3
2. minimap2=2.22
3. intervaltree=3.1.0
4. numpy=2.0.1
5. scipy=1.14.0
6. pyabpoa=1.5.2
```
## Installation
### Install by conda
```
conda
```
### Install by pip
```
pip
```
## Usage
GFHunter offers 2 steps to detect fusions: **index** and **detect**
```
GFHunter.py [-h] {index,detect} ...

positional arguments:
  {index,detect}
    index         Create index of the GFHunter
    detect        Dectect gene fusions

options:
  -h, --help      show this help message and exit
```
### Index
```
GFHunter.py index <annotationfile.gtf> <referencefile.fa> <indexdir>
```
Note: the annotation only support GENCODE currently.
### Detection
```
GFHunter.py detect <readfile.fa/fq> <indexdir>
```
| Parameter | Description | Default |
|-----:|---------------|-----|
|--middlefile|Middle file directory|./middlefile/|
|--output|Output result name|./result|
|--print_middle_output|Save the middle file after detection|NULL|
|--thread|Number of thread to use|4|
|--min_read_length|The minimum length of reads (bp) considered|50|
|--max_exon_boundary|The maximum length between breakpoint and exon boundary|30|
|--overlap_precent|The precent of overlap between reads and transcripts|0.5|
|--min_clustering_length|The minimum distance between two clusters|200|
|--least_support_reads|Least reads number to support gene fusions|2|
## Citation
## Contact
