#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 04:00:00
#SBATCH -J remove_sample_vcf
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools view and the -s flag to remove sample 101474 from the input VCF file

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/removed_101474

filename=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -s ^101474 -Oz -o "$out_dir"/"$filename"_rm.vcf.gz $1

