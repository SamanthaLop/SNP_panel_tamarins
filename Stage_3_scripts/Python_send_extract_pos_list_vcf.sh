#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J Pysend_poslist
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools pysam/0.16.0.1-python3.8.7 python/3.8.7

# This script sends the extract_pos_list_vcf to be run with python, input is the vcf file you wish to extract the positions from

scripts_dir=/proj/proj_name/nobackup/SAM/scripts
vcf_file=$1
output="$(echo $1 | cut -d . -f 1)"

python3 $scripts_dir/extract_pos_list_vcf.py $vcf_file $output
