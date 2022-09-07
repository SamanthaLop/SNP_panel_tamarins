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

# This script uses bcftools and the repeatmasker output to filter out variants in repetitive regions. Input is the interval you want "interval_2" in that 
# format, and the second input is the VCF file you wish to filter. This will only work on intervals/chromosomes 1 through 22, since those are the bed
# files that exist in the repeat_mask_dir.

module load bioinfo-tools samtools/1.12 bcftools/1.12

repeat_mask_dir=/proj/proj_name/nobackup/SAM/references/callithrix/repeatmasker/modified/

VCF=$1
interval=$(echo $(basename $1) | cut -d . -f 1 | cut -d _ -f 3,4)
repeat_mask_file="$repeat_mask_dir"/"$interval"_calJac4.fa.bed

output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/RepeatMasker_filtered

# bcftools view -T ^$repeat_dir/$REGION_calJac4.fa.out -Oz -o $output_dir/$REGION_repet_filt.vcf.gz $input_dir/$

filename_vcf=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -T ^$repeat_mask_file -Oz -o "$output_dir"/"$filename_vcf"_repet.vcf.gz $VCF

