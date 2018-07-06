# coupleNMF

Here is the source code for Integrative analysis of single cell genomics data by coupled nonnegative matrix factorizations. When different types of functional genomics data are generated on single cells from different samples of cells from the same heterogeneous population, the clustering of cells in the different samples should be coupled. We formulate this “coupled clustering” problem as an optimization problem, and propose the method of coupled nonnegative matrix factorizations (coupled NMF) for its solution. The method is illustrated by the integrative analysis of single cell RNA-seq and single cell ATAC-seq data.

For preprocessing the scATAC-seq data, please first run the preprossing script to the the openness matrix.

Example for running coupleNMF:

python coupleNMF.py -k 3 -E RAd4/E.txt -PeakO RAd4/PeakO.txt -REO RAd4/REO.txt -E_symbol RAd4/symbol.txt

For any questions about the algorithm, please contact <zduren@stanford.edu>.
