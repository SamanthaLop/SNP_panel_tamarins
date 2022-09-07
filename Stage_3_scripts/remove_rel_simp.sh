#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J remove_rel_simp
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools view and the -s flag to remove close relatives (1st degree) samples from the input vcf

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/rm_relatives
#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/rerun_combine_genotype/one_vcf_per_snp/unfiltered_rm_rels
out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/adding_to_ind_set_feb_2022/rerun_combine_genotype_gvcfs/rm_rel_simp_3

filename=$(echo $(basename $1) | cut -d . -f 1)

bcftools view -s ^101479,101480,101481 -Oz -o "$out_dir"/"$filename"_rmrelSIMP.vcf.gz $1
