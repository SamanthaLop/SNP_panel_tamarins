#!/bin/bash

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J angsd_consensus
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# Script that uses ANGSD to obtain a fasta genome from a bam file.

module load bioinfo-tools ANGSD/0.933
#module load bioinfo-tools ANGSD/0.921 this is the previous version

## Give bamfile as command input (can be run in a loop)
INPUT=$1

output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/ANGSD_consensus_seqs/101471

## Path to store raw output from angsd, will be parsed with this script so don't have to keep this
OUTPUT=./$(basename $INPUT)_angsd 

## For some reason angsd require a list of bamfiles, so I create a temporary bamlist from the command input bam
echo $INPUT > ./tmp.bamlist_$(basename $INPUT)

#running angsd with dumpcounts 3, will output allele counts for all bases at all sites where coverage > 0
#angsd -out $output_dir/"$OUTPUT" -doFasta 2 -doCounts 1 -dumpCounts 3 -b ./tmp.bamlist_$(basename $INPUT)

# running with doFasta 4 will output fasta with IUPAC codes
angsd -out $output_dir/"$OUTPUT" -doFasta 4 -doCounts 1 -dumpCounts 3 -b ./tmp.bamlist_$(basename $INPUT)
