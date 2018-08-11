#We now support mm9, mm10, hg19 and hg38 only, change the genome to your interested genome
genome=mm9  

#Download the pre-calculated files
wget http://web.stanford.edu/~zduren/CoupledNMF/Thresholding-Based%20SVD_files/common_data.tar.gz
tar -zxvf common_data.tar.gz

#Please put all the single cell ATAC-seq bed files into ./Bed/ folder
#Get PeakO.txt and PeakName.txt file
cat ./Bed/* |sort -k1,1 -k2,2n> merge.bed
ls ./Bed/ > ATAC_SampleName
macs2 callpeak -t merge.bed -f BED -n Peak --nomodel --extsize 147
cat Peak_peaks.narrowPeak|cut -f 1-3 |sort -k1,1 -k2,2n> region.bed
bedtools intersect -a region.bed -b merge.bed -wa -c -sorted>region_read_merge.bed
cat region.bed|tr '\t' '_'>region_read.bed
cat ATAC_SampleName|while read line
do
bedtools intersect -a region.bed -b ./Bed/${line} -wa -c -sorted|cut -f 4 >a
paste -d '\t' region_read.bed a>a1
cat a1 > region_read.bed
done
cat region_read.bed |cut -f 1 >PeakName.txt
n=`cat ATAC_SampleName|wc -l`
let n1=n+1
cat region_read.bed|cut -f 2-$n1 > PeakO.txt

#Get peak_gene_100k_corr.bed file
bedtools intersect -a common_data/Promoter_100k_${genome}.bed  -b region.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$3-100000-$6}'|sed 's/-//g'|sort -k1,1 -k2,2n >peak_gene_100k.bed
bedtools intersect -a peak_gene_100k.bed -b common_data/RE_gene_corr_${genome}.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{if ($4==$9) print $1,$2,$3,$4,$5,$10}'|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
