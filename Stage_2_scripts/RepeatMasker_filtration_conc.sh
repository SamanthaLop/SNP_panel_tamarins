#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J Repeat_Mask_filt
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools and the repeatmasker output to filter out variants in repetitive regions

module load bioinfo-tools samtools/1.12 bcftools/1.12

repeat_mask_file=/proj/proj_name/nobackup/SAM/references/callithrix/repeatmasker/modified/all_chr_repeat_mask.bed
output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/GATK_filtered/random_subsets_gatk_filt/concatenated
input_vcf=$1

# bcftools view -T ^$repeat_dir/$REGION_calJac4.fa.out -Oz -o $output_dir/$REGION_repet_filt.vcf.gz $input_dir/$

filename=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -T ^$repeat_mask_file -Oz -o "$output_dir"/"$filename"_repet_filt.vcf.gz $input_vcf
