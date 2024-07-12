#!/bin/bash
#conda activate allcools
chrom_size=sorted.hg38.chrom.sizes
bed_path=Dx_peaks_bed
outp=outdir

for bed_gz in `cat list.dxpeak_all`
do
bedname=${bed_gz%.sort*}

allcools generate-dataset\
	--allc_table allc_table.tsv\
	--output_path $outp\
	--chrom_size_path $chrom_size\
	--obs_dim cell \
	--cpu 1\
	--chunk_size 1\
	--regions $bedname $bed_path/$bed_gz\
	--quantifiers $bedname hyper-score CNN,CHN cutoff=0.9\
	--quantifiers $bedname hypo-score CNN,CHN cutoff=0.9
done
