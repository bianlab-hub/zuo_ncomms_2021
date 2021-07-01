#!/bin/bash

cores=8

#bedtools shuffle -i Vara_PD_CTCF_peaks.bed -g mm10.chrom.sizes > CTCF_peaks_shuffled.bed

Sertoli="/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/Sertoli_CTCF_rmdup.bw"
PreL="/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/PreL_CTCF_rmdup.bw"
LZ='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/LZ_CTCF_rmdup.bw'
PD='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/PD_CTCF_rmdup.bw'
Vara_PD='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/GSM3840086_CovPDCTCF.bw'

CTCF_only='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/Vara_PD_CTCF_peaks.bed'
CTCF_REC8_both='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/Vara_PD_CTCF_peaks_Vara_PD_REC8_peaks.bed'
REC8_only='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/Vara_PD_REC8_peaks.bed'
control='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/CTCF_peaks_shuffled.bed'



computeMatrix scale-regions -S ${Sertoli} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${control} -o Sertoli.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m Sertoli.mat.gz -out Sertoli_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m Sertoli.mat.gz -out Sertoli_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

computeMatrix scale-regions -S ${PreL} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${control} -o PreL.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m PreL.mat.gz -out PreL_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m PreL.mat.gz -out PreL_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

computeMatrix scale-regions -S ${LZ} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${control} -o LZ.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m LZ.mat.gz -out LZ_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m LZ.mat.gz -out LZ_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

computeMatrix scale-regions -S ${PD} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${control} -o PD.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m PD.mat.gz -out PD_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m PD.mat.gz -out PD_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

computeMatrix scale-regions -S ${Vara_PD} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${control} -o Vara_PD.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m Vara_PD.mat.gz -out Vara_PD_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m Vara_PD.mat.gz -out Vara_PD_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

#bedtools shuffle -i Vara_PD_REC8_peaks.bed -g mm10.chrom.sizes > REC8_peaks_shuffled.bed
Vara_PD_REC8='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/GSM3840088_CovPDREC8.bw'
REC8_control='/bianlab/zuowu/meiosis/Cut/alignment_trim/bigwig/deeptools/REC8_peaks_shuffled.bed'


computeMatrix scale-regions -S ${Vara_PD_REC8} -R ${CTCF_only} ${CTCF_REC8_both} ${REC8_only} ${REC8_control} -o Vara_PD_REC8.mat.gz -b 2000 -a 2000 -m 2000 --skipZeros -p $cores &
plotHeatmap -m Vara_PD_REC8.mat.gz -out Vara_PD_REC8_heatmap.pdf --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotProfile -m Vara_PD_REC8.mat.gz -out Vara_PD_REC8_profile.pdf --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --plotHeight 10 --plotWidth 12 &

plotHeatmap -m Vara_PD_REC8.mat.gz -out Vara_PD_REC8_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m Sertoli.mat.gz -out Sertoli_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m PreL.mat.gz -out PreL_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m LZ.mat.gz -out LZ_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m PD.mat.gz -out PD_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m Vara_PD.mat.gz -out Vara_PD_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
plotHeatmap -m Vara_PD_REC8.mat.gz -out Vara_PD_REC8_heatmap.png --sortUsing sum --dpi 300 -T "" --regionsLabel "CTCF_only" "CTCF_REC8_both" "REC8_only" "control" --startLabel ""  --endLabel "" --legendLocation "upper-right" --heatmapWidth 12 &
