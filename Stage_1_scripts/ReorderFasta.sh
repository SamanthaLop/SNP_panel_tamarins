#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 02:00:00
#SBATCH -J ReorderFasta
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=ALL

# this script reorders the raw sequencing files so that the paired end reads match to one another. It takes two arguments, *_1.fastq.gz and *_2.fastq.gz
# It will generate three output files: *_1_sorted.fastq.gz, *_2_sorted.fastq.gz and *_singletons.fastq.gz

module load bioinfo-tools bbmap

# write this in the terminal
# sbatch ReorderFasta.sh /proj/proj_name/nobackup/TAMARINS/RAW_READS/101455_1.fastq.gz /proj/proj_name/nobackup/TAMARINS/RAW_READS/101455_2.fastq.gz

file1=$1
file2=$2
OUTPUT=/proj/proj_name/nobackup/SAM/reordered_tamarins

name=$(basename $file1)
SM=$(echo $name | cut -d _ -f 1)

repair.sh in1=$file1 in2=$file2 out1=$OUTPUT/"$SM"_1_sorted.fastq.gz out2=$OUTPUT/"$SM"_2_sorted.fastq.gz outs=$OUTPUT/"$SM"_singletons.fastq.gz repair
