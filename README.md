# coupleNMF

## Introduction
Here is the source code for Integrative analysis of single cell genomics data by coupled nonnegative matrix factorizations. When different types of functional genomics data are generated on single cells from different samples of cells from the same heterogeneous population, the clustering of cells in the different samples should be coupled. We formulate this “coupled clustering” problem as an optimization problem, and propose the method of coupled nonnegative matrix factorizations (coupled NMF) for its solution. The method is illustrated by the integrative analysis of single cell RNA-seq and single cell ATAC-seq data.

## Preprocessing
For preprocessing the scRNA-seq data, please following the standard processing pipline to get the expression matrix, where each row represents a gene, each column represents a cell.

For preprocessing the scATAC-seq data, please first put all the .bam files for each cell into a folder. Then run the preprossing script we provided to get the the openness matrix, REO and PeakO. 

## Running coupleNMF
coupleNMF receives 7 parameters:

-k         the clustering numbers

-PeakO     the location of PeakO matrix

-REO       the location of REO matrix

-E         the location of E matrix

-E_symbol  the location of gene symbol file

-s         species (human or mouse)

-ref       the reference genome (mm9, mm10, hg19, hg38)  


### Example

```
python coupleNMF.py -k 3 -E RAd4/E.txt -PeakO RAd4/PeakO.txt -REO RAd4/REO.txt -E_symbol RAd4/symbol.txt -s mouse -ref mm9

```


## Requirements
* sklearn
* pandas
* scipy
* itertools
* argparse 
* itertools
* MACS

For any questions about the algorithm, please contact <zduren@stanford.edu>.
