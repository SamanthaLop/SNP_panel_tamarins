#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -J 2bcf_filter
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools to further filter the variants.

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/bcf2_filtered

filename=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -i 'INFO/MAF>0.1' --min-alleles 2 --max-alleles 2 --types "snps"  -f 'PASS' -Oz -o "$out_dir"/"$filename"_bcf2.vcf.gz $1

#bcftools view -i 'INFO/MAF>0.2' --min-alleles 2 --max-alleles 2 --types "snps" -f 'PASS' -Oz -o "$out_dir"/"$filename"_bcf2maf02.vcf.gz $1
