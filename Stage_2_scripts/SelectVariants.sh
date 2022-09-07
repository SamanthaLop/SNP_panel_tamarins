#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J SelectVariants
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# script uses GATK to remove variants that have at least one sample with genotype set to no call (uses the AB filtered output)

module load bioinfo-tools GATK/4

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls

filename=$(echo $(basename $1) | cut -d . -f 1)

# gatk IndexFeatureFile -I $1
gatk SelectVariants -R $ref_genome -V $1 --max-nocall-number 0 -O "$output_dir"/"$filename"_rmnocall.vcf.gz
