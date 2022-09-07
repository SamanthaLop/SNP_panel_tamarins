#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J VCFtools_stats
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

## script to obtain vcf stats from a subsetted vcf file.
module load bioinfo-tools bcftools vcflib vcftools

SUBSET_VCF=$1
stats_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/GATK_filtered/random_subsets_gatk_filt/concatenated
OUT=$(echo $(basename $1) | cut -d . -f 1)

# To calculate Allele frequency
vcftools --gzvcf $SUBSET_VCF --freq2 --out "$stats_dir"/$OUT --max-alleles 2
# output file ends with .frq

# To calculate mean depth of coverage per individual
vcftools --gzvcf $SUBSET_VCF --depth --out "$stats_dir"/$OUT
# output file ends with .idepth

# To calculate mean depth of coverage for each site
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out "$stats_dir"/$OUT
# output file ends with .ldepth-mean

# To calculate site quality score for each site
vcftools --gzvcf $SUBSET_VCF --site-quality --out "$stats_dir"/$OUT
# output file ends with .lqual

# To calculate the proportion of missing data per sample
vcftools --gzvcf $SUBSET_VCF --missing-indv --out "$stats_dir"/$OUT
# output file ends with .imiss

# To  calculate the proportion of missing data per site
vcftools --gzvcf $SUBSET_VCF --missing-site --out "$stats_dir"/$OUT
# output file ends with .lmiss

# To calculate heterozygosity and inbreeding coefficient (F) per individual
vcftools --gzvcf $SUBSET_VCF --het --out "$stats_dir"/$OUT
# output file ends with .het

# extra, to calculate the depth for each genotype. This will be a very large file
vcftools --gzvcf $SUBSET_VCF --geno-depth --out "$stats_dir"/$OUT
# output file ends with .gdepth

# extra, to calculate depth per site summed across all individuals!
vcftools --gzvcf $SUBSET_VCF --site-depth --out "$stats_dir"/$OUT
# output file ends with .ldepth

# to calculate p value for each site from a HWE, also obsv and exp het
vcftools --gzvcf $SUBSET_VCF --hardy --out "$stats_dir"/$OUT
