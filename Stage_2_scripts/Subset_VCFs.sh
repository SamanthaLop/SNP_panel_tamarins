#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J Subset_VCFs
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This file takes the vcf file (resulting from GenotypeGVCFs) and subsets a given number of variants from it to be able to explore different
# filtering methods or concatenate the subsets to apply the kinship analysis with KING.

module load bioinfo-tools bcftools vcflib GATK

#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/subsets
# dir for snakemake pipeline troubleshooting:
out_dir=/proj/proj_name/nobackup/SAM/snakemake_attempt/genotypeGVCFs/subsets


#calculate r needed to obtain either x (100,000) or y (10,000) variants

#total=$(echo $(bcftools view -H $1 | wc -l))
total_snps=$(echo $(bcftools view -H $1 --types 'snps' | wc -l))
echo $total_snps

#x=10000
# for snakemake pipeline troubleshooting:
#x=1000
x=100

r=$(echo "scale=6; $x/$total_snps" | bc -l)

file=$(basename $1)
file_name=$(echo $file | cut -d . -f 1)

# make the subset with the calulated r and subsetting only variant sites (SNPs) NOTE: this script will fail if you already have a file under that name
# so be sure to get rid of your previous files before running
#bcftools view $1 --types 'snps' | vcfrandomsample -r $r >"$out_dir"/"$file_name"_subset.vcf
bcftools view $1 --types 'snps' | vcfrandomsample -r $r >"$out_dir"/"$file_name".vcf

# compress subset
#bgzip "$out_dir"/"$file_name"_subset.vcf
bgzip "$out_dir"/"$file_name".vcf


# index subset
#bcftools index "$out_dir"/"$file_name"_subset.vcf.gz

gatk IndexFeatureFile -I "$out_dir"/"$file_name".vcf.gz
