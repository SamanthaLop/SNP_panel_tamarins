#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J Pysend_compare
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools pysam/0.16.0.1-python3.8.7 python/3.8.7

# This script sends the compare_positions_vcf.py script to python to be run. Input should simply be the interval you want it to run on.
# This can be run like this: sbatch Python_send_compare_positions_vcf.sh interval_2 (including the underscore)
# The script will call the simp position list and the lwed position list generated with extract_pos_list_vcf.py

scripts_dir=/proj/proj_name/nobackup/SAM/scripts

REGION=$1

# 0.30-0.45 MAF female het snps
simp_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/x_chromosome/NEW_MAF
lwed_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/x_chromosome/NEW_MAF

simp_pos_list="$simp_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmlwedd_rmrelSIMP_rmmalesSIMP_vcffixup_MAFfilt5_pos_list.tsv
lwed_pos_list="$lwed_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmsimp_rmrelLWED_rmmalesLWED_vcffixup_MAFfilt5_pos_list.tsv

output=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/x_compared/NEW_MAF_COMPARED/"$REGION"_MAF0.30_0.45

python3 $scripts_dir/compare_positions_vcf.py $simp_pos_list $lwed_pos_list $output


#  0.30-0.45 MAF individual identification
#simp_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/MAFfilt_updated
#lwed_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/MAFfilt_updated

#simp_pos_list="$simp_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmlwedd_rmrelSIMP_vcffixup_MAFfilt5_pos_list.tsv
#lwed_pos_list="$lwed_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmsimp_rmrelLWED_vcffixup_MAFfilt5_pos_list.tsv

#output=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/shared_variants/MAF0.30_0.45/"$REGION"_MAF0.30_0.45

#python3 $scripts_dir/compare_positions_vcf.py $simp_pos_list $lwed_pos_list $output


# 0.15-0.35 MAF
#simp_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/MAFfilt/poslist/
#lwed_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/MAFfilt/poslist/

#simp_pos_list="$simp_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmlwedd_rmrelSIMP_vcffixup_MAFfilt2_pos_list.tsv
#lwed_pos_list="$lwed_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmsimp_rmrelLWED_vcffixup_MAFfilt2_pos_list.tsv

#output=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/shared_variants/"$REGION"_2

#python3 $scripts_dir/compare_positions_vcf.py $simp_pos_list $lwed_pos_list $output

# ---------------x chromosome
#simp_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/x_chromosome
#lwed_pos_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/lwedd/x_chromosome

#simp_pos_list="$simp_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmlwedd_rmrelSIMP_rmmalesSIMP_vcffixup_MAFfilt3_pos_list.tsv
#lwed_pos_list="$lwed_pos_dir"/gatk_filtered_"$REGION"_rm_repet_bcf1_bcf2_AB_rmnocall_rmsimp_rmrelLWED_rmmalesLWED_vcffixup_MAFfilt3_pos_list.tsv

#output=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/x_compared/"$REGION"_4

#python3 $scripts_dir/compare_positions_vcf.py $simp_pos_list $lwed_pos_list $output
