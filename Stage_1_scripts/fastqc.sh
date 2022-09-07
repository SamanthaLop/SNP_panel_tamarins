#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 05:00:00
#SBATCH -J fastqc
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=ALL

# This script runs fastqc analysis on raw reads (input is fastq)

# Load the modules
module load bioinfo-tools FastQC

# -t 20 means it can multi-thread, it uses several cores at a time to be more efficient.
# -o specifies the output

fastqc /proj/proj_name/nobackup/TAMARINS/sorted_tamarins/$1 -t 20 -o /proj/proj_name/nobackup/SAM/fastqc/
