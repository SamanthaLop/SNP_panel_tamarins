#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J MAF_filtering
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools to filter VCF files to exclude variants with MAF >0.45 and <0.30
# remember, if you are running this on a vcf where you have removed samples, you need to run the vcf through vcffix
# to update INFO fields!!

module load bioinfo-tools bcftools

#REMEMBER TO UPDATE OUTDIR!

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/x_chromosome/NEW_MAF

filename=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -i 'INFO/MAF>0.30 && INFO/MAF<0.45'  -f 'PASS' -Oz -o "$out_dir"/"$filename"_MAFfilt5.vcf.gz $1
