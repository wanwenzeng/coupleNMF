#Please put all the single cell RNA-seq bed files into ./Bed/ folder

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
cat region_read.bed |cut -f 1 >region_name.txt
cat RE_v2.bed|tr '\t' '_'>RE_read.bed
cat ATAC_SampleName|while read line
do
bedtools intersect -a RE_v2.bed -b ./Bed/${line} -wa -c|cut -f 4 >a
paste -d '\t' RE_read.bed a>a1
cat a1 > REO.txt
done
