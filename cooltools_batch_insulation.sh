
#!/bin/bash
conda activate base

threads=12


for filename in *_10kb.cool; do 
	sample=$(basename "$filename" .cool)
	echo $sample
	cooler balance ${sample}.cool -p $threads
	cooltools compute-expected ${sample}.cool -o ${sample}_cis_expected.tsv -t cis -p $threads 
	cooltools compute-expected ${sample}.cool -o ${sample}_trans_expected.tsv -t trans -p $threads
	cooltools diamond-insulation ${sample}.cool 500000	> ${sample}_insulation_500kb.tsv
done