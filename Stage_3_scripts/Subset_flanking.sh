#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:59:00
#SBATCH -J Subset_flanking
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools vcflib

# This script was written to be able to subset a list of potential SNPs (obtained from fixeddifferences) from a vcf so that I can use the --thin flag from vcftools
# I need to feed this tool a vcf file, so I get the large filtered VCF file  for the chromosome and extract the potential SNPs. I then apply the --thin option 
# to this vcf and together with bcftools and awk, I get a list of positions that are at least n bases apart from each other.

REGION=$1
input_vcf=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz
flanking_sites=/proj/proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_speciesSNP_pos_list.tsv
filename=$(echo $(basename $input_vcf) | cut -d . -f 1)
output_dir=/proj/proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/vcfs


bcftools view --regions-file $flanking_sites -Oz --output "$output_dir"/"$filename"_VCF.vcf.gz $input_vcf
