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

**coupleNMF outputs 3 files:**

* scATAC-result.txt                       the clustering results for scATAC-seq

* scRNA-result.txt                        the clustering results for scRNA-seq

* cluster-specific-peaks-genes-pairs.txt  the cluster-specific peak-gene pairs. First column is the gene name, second column is the peak name, the third column is the p-value for gene and last column is the p-value for peak. 



### Example

```
python coupleNMF.py -k 2 -E exampledata/E.txt -PeakO exampledata/PeakO.txt -E_symbol exampledata/symbol.txt -P_symbol exampledata/PeakName.txt -pe exampledata/peak_gene_100k_corr.bed  -lambda1 0.04 -lambda2 25

```


## Connection of Cell Ranger

Since the Cell Ranger （https://github.com/10XGenomics/cellranger) is one of the most popular single cell preprocessing software, our algorithm can further process the output of Cell Ranger. For example, one of the Cell Ranger output is the matrix.mtx as

```
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"format_version": 2, "software_version": "1.2.0"}
579 668 382948
1 1 6622
2 1 14970
3 1 58
4 1 339
5 1 749
6 1 1014
…………
…………
```

Then we can use the following easy script to transfer the output matrix.mtx file to the input format file our algorithm

```
from scipy.io import mmread
import numpy as np

mat = mmread("matrix.mtx")
np.savetxt('matrix.txt', mat.todense(),  delimiter='\t', fmt="%d")
```


## Requirements
* sklearn
* pandas
* scipy
* itertools
* argparse 
* itertools
* MACS
* umap (https://github.com/lmcinnes/umap)



**For any questions about running the software, please contact <zengww14@mails.tsinghua.edu.cn>.**

**For any questions about the algorithm, please contact <zduren@stanford.edu>.**
