#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:29:00
#SBATCH -J vcftools_thin
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools bcftools vcflib vcftools

# This script has the commands I used to be able to extract a list of SNPs that are sufficiently far apart from each other (500,000bp) using  the --thin flag from vcftools
# The script couldn't run for some reason so I just copied the commands I used for posterity.

input_vcf=$1
out_dir=$(echo $(realpath $1))
filename=$(echo $(basename $1) | cut -d . -f 1)
REGION=$(echo $(basename $1) | cut -d . -f 3,4)

# ran this on the terminal, not as script. Kept failing for some reason.
#vcftools --gzvcf $input_vcf --thin 500000 --recode --out "$output_dir"/"$filename"_thin.vcf.gz
vcftools --gzvcf $input_vcf --thin 100000 --recode --out "$output_dir"/"$filename"_thin.vcf.gz

# ran this on the terminal, not as a script. Kept failing for some reason.
bcftools view -H "$output_dir"/"$filename"_thin.vcf.gz | awk '{print $1"\t"$2}' > "$output_dir"/"$REGION"_thin_pos_list.tsv
