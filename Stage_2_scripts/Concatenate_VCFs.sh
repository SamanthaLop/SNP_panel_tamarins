#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J Concatenate_VCFs
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools concatenate to merge or combine vcf files to get a single output file. I will use it to concatenate the  subsets
# from each chromosome (autosomal). Make sure the filenames file contains all the filenames you wish to merge.

module load bioinfo-tools bcftools vcflib

list=/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/candidate_snps_vcfs/filenames.txt
out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/candidate_snps_vcfs

cd $out_dir

bcftools concat --file-list $list --naive -O z > "$out_dir"/concatenated_candidate_snps_indID.vcf.gz

# -O z outputs compressed VCF file
# --naive concatenates VCF or BCF files without recompression. It's very fast but requires all files to be the same type (VCF or BCF) and have the same headers.
