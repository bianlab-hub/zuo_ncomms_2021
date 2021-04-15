#For each bed file containing locations of TAD boundaries, use bedtools slop to obtain the 1Mb regions centered at the boundaries 
#Use coolpup.py --local to perform the pileup, the expected files used for normalization were generataed using cooltools
#coolpup.py /data05/combined_spermatocytes/sertoli/fragment_db/sertoli_combined_10kb.cool PD_CTCF_REC8_overlapping_10kb.bed --nshifts 10 --pad 500 --mindist 1000000 --maxdist 2000000 --outname sertoli_loop_1M_2M.txt --log WARNING &
#--expected ../expected/${stage}_10kb_cis_expected.tsv 

for stage in I Sertoli PreL L Z P D
do
echo $stage
coolpup.py ../${stage}_10kb.cool PD_CTCF_REC8_overlapping_10kb.bed --nshifts 10 --pad 500 --mindist 1000000 --maxdist 2000000 --pad 500 --outname ${stage}_CTCF_REC8_pileup_1m_2m.txt --n_proc 8 --log WARNING
coolpup.py ../${stage}_10kb.cool PD_CTCF_REC8_overlapping_10kb.bed --expected ../expected/${stage}_10kb_cis_expected.tsv --pad 500 --mindist 1000000 --maxdist 2000000 --pad 500 --outname ${stage}_CTCF_REC8_pileup_1m_2m_expected.txt --n_proc 8 --log WARNING 
coolpup.py ../${stage}_10kb.cool PD_CTCF_REC8_overlapping_10kb.bed --nshifts 10 --pad 500 --mindist 500000 --maxdist 2000000 --pad 500 --outname ${stage}_CTCF_REC8_pileup_500k_2m.txt --n_proc 8 --log WARNING
coolpup.py ../${stage}_10kb.cool PD_CTCF_REC8_overlapping_10kb.bed --expected ../expected/${stage}_10kb_cis_expected.tsv --pad 500 --mindist 500000 --maxdist 2000000 --pad 500 --outname ${stage}_CTCF_REC8_pileup_500k_2m_expected.txt --n_proc 8 --log WARNING 
 
done



