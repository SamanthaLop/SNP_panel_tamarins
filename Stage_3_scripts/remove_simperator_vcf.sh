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

# This script uses bcftools view and the -s flag to remove Saguinus imperator samples from the input vcf

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/candidate_snps_vcfs_150/test_hwe/lwed
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/adding_to_ind_set_feb_2022/HWE/only_snps_vcfs
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/redoing_sex_snps_per_species

filename=$(echo $(basename $1) | cut -d . -f 1)

#bcftools view -s ^101467,101468,101477,101478,101479,101480,101481 -Oz -o "$out_dir"/"$filename"_rmsimp.vcf.gz $1

# below is what needs to be removed if you've already removed related individuals (if you try to remove a sample not present in the vcf, bcftools
# will throw and error. This is why I did this. 
bcftools view -s ^101467,101468,101477,101478 -Oz -o "$out_dir"/"$filename"_rmsimp.vcf.gz $1

