# Gene Expression Data Files
This directory should contain the gene expression data files produced from the Cellranger pipeline, divided into separate folders based on sample. Each folder fhould contain the three following files: barcodes.tsv, genes.tsv, matrix.mtx. 

The organization of this directory should look like this: 
```bash
gene_expression
├── README.md
├── RA1
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── RA2
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
├── RA3
│   ├── barcodes.tsv
│   ├── genes.tsv
│   └── matrix.mtx
└── RA4
    ├── barcodes.tsv
    ├── genes.tsv
    └── matrix.mtx
```
