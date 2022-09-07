#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J Python_send
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# this script is to send the AB_Filter.py as a job to uppmax

module load bioinfo-tools pysam/0.16.0.1-python3.8.7

REGION=$(echo $(basename $1) | cut -d . -f 1)

scripts_dir=/proj/proj_name/nobackup/SAM/scripts
bcf2_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/bcf2_filtered
output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered

python3 "$scripts_dir"/AB_Filter.py "$bcf2_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2.vcf.gz "$output_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB.vcf.gz


