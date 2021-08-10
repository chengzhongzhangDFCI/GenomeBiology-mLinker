#!/bin/bash

declare -a chromosomes=('chr21')

## FOR NA12878_WG_210 ########
##############################
suffix='NA12878_WG_210_0317_test'
vcf='/czlab/Data/NA12878/10X/HaplotypeCalling/NA12878.hets.chr1-X.BA_withRef.vcf.gz'
#vcf='/czlab/Data/NA12878/10X/HaplotypeCalling/NA12878.hets.chr1-X.hets.vcf.gz'
bam='/czlab/Data/NA12878/10X/NA12878_WGS_210.bam'
technology='tenx'
kmersize=20
############################################


refgendir='/czlab/References/GRCh38/hg38_chroms_fa/'

for chr in ${chromosomes[@]}
do
	nohupname='nohup_extract_'${suffix}'_'${chr}'.out'
	nohup time ./mlinker extract -v ${vcf} -i ${bam} -c ${chr} -e ${technology} -d ${refgendir} -l 100 -k ${kmersize} -n ${suffix}  &> ${nohupname} &
done
