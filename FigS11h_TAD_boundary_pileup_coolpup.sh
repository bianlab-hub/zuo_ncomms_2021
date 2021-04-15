#For each bed file containing locations of TAD boundaries, use bedtools slop to obtain the 1Mb regions centered at the boundaries 
#Use coolpup.py --local to perform the pileup, the expected files used for normalization were generataed using cooltools

for filename in *_10kb_TAD_boundary_0.1.bed; do 
	sample=$(basename "$filename" _10kb_TAD_boundary_0.1.bed)
	echo $sample
	bedtools slop -i ${sample}_10kb_TAD_boundary_0.1.bed -g mm10.chrom.sizes -b 500000 > ${sample}_10kb_TAD_boundary_plus_minus_500kb.bed
	coolpup.py ../${sample}_10kb.cool ${sample}_10kb_TAD_boundary_plus_minus_500kb.bed --expected ../expected/${sample}_10kb_cis_expected.tsv --local --pad 500 --outname ${sample}_TAD_boundary_pileup.txt --n_proc 8 --log WARNING 
done



