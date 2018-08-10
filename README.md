# coupleNMF

## Introduction
Here is the source code for Integrative analysis of single cell genomics data by coupled nonnegative matrix factorizations. When different types of functional genomics data are generated on single cells from different samples of cells from the same heterogeneous population, the clustering of cells in the different samples should be coupled. We formulate this “coupled clustering” problem as an optimization problem, and propose the method of coupled nonnegative matrix factorizations (coupled NMF) for its solution. The method is illustrated by the integrative analysis of single cell RNA-seq and single cell ATAC-seq data.

## Preprocessing
For preprocessing the scRNA-seq data, please following the standard processing pipline to get the expression matrix, where each row represents a gene, each column represents a cell.

For preprocessing the scATAC-seq data, please first put all the .bam files for each cell into a folder. Then run the preprossing script we provided to get the the openness matrix, PeakO and PeakName. 


## Running coupleNMF
**coupleNMF receives 8 parameters:**

* -k         the clustering numbers

* -PeakO     the location of PeakO matrix

* -E         the location of E matrix

* -E_symbol  the location of gene symbol file

* -P_symbol  the location of peak symbol file

* -pe        the location of pre-calculated peak-gene interactions file 

* -lambda1   the hyper-paramters lambda1 to control term of the NMF for E 

* -lamdba2   the hyper-paramters lambda2 to control coupled term

Note:-k, -PeakO, -E, -E_symbol, -P_symbol, -pe are the must-have parameters; 
-lambda1, -lamdba2 are optional parameters. If coupleNMF does not receive -lambda1 and -lambda2, it will choose the best parameters automatically.



### Example

```
python coupleNMF.py -k 2 -E exampledata/E.txt -PeakO exampledata/PeakO.txt -E_symbol exampledata/symbol.txt -P_symbol exampledata/PeakName.txt -pe common_data/peak_gene_100k_corr.bed  -lambda1 25 -lambda2 0.04

```


## Requirements
* sklearn
* pandas
* scipy
* itertools
* argparse 
* itertools
* MACS




**For any questions about running the software, please contact <zengww14@mails.tsinghua.edu.cn>.**

**For any questions about the algorithm, please contact <zduren@stanford.edu>.**
