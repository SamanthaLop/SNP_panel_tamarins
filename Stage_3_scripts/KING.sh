#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J KING
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script will male a PLINK binary format set of files from a VCF input file and perform KING kinship analysis.

module load bioinfo-tools plink

# REMEMBER TO UPDATE THE OUTPUT DIRECTORY AS NEEDED
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/subsets/KING/simperator_king
subset_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/subsets/KING
file=$(basename $1)
file_name=$(echo $file | cut -d . -f 1)

plink --vcf $file --make-bed --out "$subset_dir"/"$file_name"_plink --allow-extra-chr --maf 0.05

# I get an error if I don't add the "--allow-extra-chr" option.

# replace chromosome IDs with integers in the bim file so that KING can run. I'm sure there's a smarter way to do this but this will work for now.

sed -i 's/CM018917.1/1/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018918.1/2/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018919.1/3/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018920.1/4/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018921.1/5/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018922.1/6/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018923.1/7/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018924.1/8/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018925.1/9/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018926.1/10/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018927.1/11/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018928.1/12/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018929.1/13/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018930.1/14/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018931.1/15/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018932.1/16/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018933.1/17/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018934.1/18/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018935.1/19/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018936.1/20/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018937.1/21/g' "$subset_dir"/"$file_name"_plink.bim
sed -i 's/CM018938.1/22/g' "$subset_dir"/"$file_name"_plink.bim

/home/samlope/glob/king -b "$subset_dir"/"$file_name"_plink.bed --kinship 


