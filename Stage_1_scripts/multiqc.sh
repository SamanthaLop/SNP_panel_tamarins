#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 05:00:00
#SBATCH -J multiqc
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=ALL

# This script uses the output of fastqc to make a report summary of all samples. It needs no input, just makes sure the correct directory (the one with the fastqc
# output) is indicated below.

module load bioinfo-tools
module load MultiQC

#multiqc /proj/proj_name/nobackup/SAM/fastqc/*.zip
multiqc /proj/proj_name/nobackup/SAM/qualimap/
