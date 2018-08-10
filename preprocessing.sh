#Please put all the single cell RNA-seq bed files into ./Bed/ folder

#Get PeakO.txt file
cat ./Bed/* |sort -k1,1 -k2,2n> merge.bed
ls ./Bed/ > ATAC_SampleName

module load gcc
ml load python
macs2 callpeak -t merge.bed -f BED -n Peak --nomodel --extsize 147
cat Peak_peaks.narrowPeak|cut -f 1-3 |sort -k1,1 -k2,2n> region.bed
bedtools intersect -a region.bed -b merge.bed -wa -c -sorted>region_read_merge.bed
cat region.bed|tr '\t' '_'>region_read.bed
cat ATAC_SampleName|while read line
do
bedtools intersect -a region.bed -b ./Bed/${line} -wa -c -sorted|cut -f 4 >a
paste -d '\t' region_read.bed a>a1
cat a1 > PeakO.txt
done
#Get peak_gene_100k_corr.bed file
genome=mm9
bedtools intersect -a common_data/Promoter_100k_${genome}.bed  -b Peak.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$3-100000-$6}'|sed 's/-//g'|sort -k1,1 -k2,2n >peak_gene_100k.bed
bedtools intersect -a peak_gene_100k.bed -b common_data/RE_gene_corr_${genome}.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{if ($4==$9) print $1,$2,$3,$4,$5,$10}'|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
