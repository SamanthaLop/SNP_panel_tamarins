#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J remove_males
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools view and the -s flag to remove male samples from the input vcf

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/x_chromosome
#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/x_chromosome

filename=$(echo $(basename $1) | cut -d . -f 1)

#remove leontocebus males
#bcftools view -s ^101455,101459,101466,101470,101473 -Oz -o "$out_dir"/"$filename"_rmmalesLWED.vcf.gz $1

# remove simperator males
bcftools view -s ^101467 -Oz -o "$out_dir"/"$filename"_rmmalesSIMP.vcf.gz $1
