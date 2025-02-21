# How to use the evaluation
Just keep the "evaluation/" structure, we can use the code in "evaluation/code/" to evaluate the result of GFHunter and four other tools.
```
evaluation
├─code
└─result
    ├─All-Simulation
    │  ├─Nagetive
    │  ├─ONT10x
    │  ├─ONT20x
    │  ├─ONT30x
    │  ├─ONT50x
    │  ├─PB10x
    │  ├─PB20x
    │  ├─PB30x
    │  └─PB50x
    ├─HCT-116
    ├─MCF-7
    ├─NGS
    ├─Non-tumor-data
    ├─Output
    │  ├─Fig.2 e&f
    │  ├─Fig.2 g
    │  ├─Fig.3 c1
    │  ├─Fig.3 c2
    │  ├─Fig.3 d
    │  ├─Fig.3 e
    │  ├─Fig.4 a
    │  ├─Fig.4 b
    │  └─Table 1
    ├─SKBR-3
    └─gencode.v47.annotation.gtf # not in github
```
- Note: user should download an annotation: gencode.v47.annotation.gtf from genecode as follow.
```
cd result
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```
