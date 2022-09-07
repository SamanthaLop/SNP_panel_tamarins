#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J Pysend_fixedd
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools pysam/0.16.0.1-python3.8.7 python/3.8.7

# FixedDifferences.py should be run like this:
# python3 fixedDifferences.py -i /your/vcf.vcf -o output.tsv --popfile sample_species.txt

scripts_dir=/proj/proj_name/nobackup/SAM/scripts

REGION=$(echo $(basename $1) | cut -d . -f 1)

vcf_file=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz
output=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/fixeddifferences/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall.tsv
popfile=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/fixeddifferences/sample_species.txt

python3 $scripts_dir/fixeddifferences.py -i $vcf_file -o $output --popfile $popfile
