#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 6-00:00:00
#SBATCH -J HaplotypeCaller
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script runs haplotypecaller. Haplotypecaller should be run on every combination of sample x chromosome. It is useful to invoke this
# this script with a loop such as this:

# for i in /proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/*.bam
#       for y in /proj/proj_name/nobackup/SAM/references/callithrix/CHROMOSOMES/*
#       do sbatch HaplotypeCaller.sh $i $y
#       done

module load bioinfo-tools samtools picard GATK/4

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna

out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/HAPLOTYPECALLER

name1=$(basename $1)
SM=$(echo $name1 | cut -d _ -f 1| cut -d . -f 2)
interval=$(basename $2)
interval_number=$(echo $interval | cut -d . -f 1)

gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R $ref_genome \
	-I $1 \
	-L $2 \
	-O "$out_dir"/"$SM"_"$interval_number".g.vcf.gz \
	-ERC GVCF


