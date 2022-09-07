#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J obtain_flanking
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools BEDTools/2.29.2

# script uses bedtools to extract n (-b flag) position ranges before and after SNP. Input is a bed file with three columns: chromosome name (not chr1, the
# name has to be the one used in the genome, contains letters and numbers, can be found in the NCBI website for the genome), position, and in the third column
# we can have either the position next to the SNP (position + 1) or just the same position as the SNP, such that the second and third columns are the same. If
# the second and third columns are the same, then the output ranges will NOT include the SNP. Since we aim to filter the flanking regions now and not the SNP,

genome_file=/proj/proj_name/nobackup/SAM/references/callithrix/cal_jac_genome_file.txt
input_list_positions=$1

output=$(echo $(basename $1) | cut -d . -f 1)

# update this when running again!
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/lists/350_flanking_lists
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/individual_species_specific_SNPs/lwed/lists/1000_flanking_intervals
out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/lists

bedtools flank -i $input_list_positions -g $genome_file -b 400 > "$out_dir"/"$output"_flanking400.bed
