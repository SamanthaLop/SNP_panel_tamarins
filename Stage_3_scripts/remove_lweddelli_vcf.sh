#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J remove_sample_vcf
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools view and the -s flag to remove Leontocebus weddelli samples from the input vcf

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/

filename=$(echo $(basename $1) | cut -d . -f 1)
#filename=$(echo $(basename $1) | cut -d . -f 1,2,3,4)

#bcftools view -s ^101455,101456,101457,101458,101459,101460,101461,101462,101463,101465,101466,101469,101470,101471,101472,101473,101475,101476 -Oz -o "$out_dir"/"$filename"_rmlwedd.vcf.gz $1

#in case you've already remove relative individuals, use this:
bcftools view -s ^101455,101458,101459,101460,101461,101465,101466,101470,101472,101473 -Oz -o "$out_dir"/"$filename"_rmlwedd.vcf.gz $1
