#!/bin/bash -l

#SBATCH -A PROJ_DETAILS
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -J bwa_tamarin
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# Adapted from a script written by Axel Jensen
# This script follows RG_data_script1.sh, maps the reads to the reference genome and produces a coordinate-sorted bam file

module load bioinfo-tools bwa samtools GATK/4 picard

# output directory
outdir=/proj/proj_name/nobackup/SAM/tamarins_mapped

# reference genome
ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna

# raw reads directory
reads_dir=/proj/proj_name/nobackup/TAMARINS/sorted_tamarins

# scripts directory
scripts=/proj/proj_name/nobackup/SAM/scripts

name=$(basename $1)
SM=$(echo $name | cut -d _ -f 1)

# PU = Platform unit (flowcell barcode, lane, sample barcode)
PU=$(zcat "$reads_dir"/"$SM"_1_sorted.fastq.gz | head -n 1 | cut -d ":" -f 2-4 | cut -d "@" -f 2 | sed "s|":"|"."|g")

# the next line	sends the next script to be run	AFTER this script is finished.
sbatch -J $(echo "$SM".merge_script3 | cut -d "_" -f 1-2) --dependency=afterok:"$SLURM_JOB_ID" "$scripts"/merge_bams_script3.sh "$outdir"/mi."$SM"_"$PU".u.bam "$outdir"/sorted."$SM"_"$PU".bam

# ----- STEP 4 - map and pipe to sortsam to get coordinate-sorted bamfile as output
bwa mem -M -t 20 -p $ref_genome $1 | samtools view -Sb -q 30 - | \
java -Xmx6G -jar "$PICARD_HOME"/picard.jar SortSam \
        -INPUT /dev/stdin \
        -OUTPUT "$outdir"/sorted."$SM"_"$PU".bam \
        -SORT_ORDER coordinate \
        -CREATE_INDEX true

rm -f "$outdir"/"$SM"_"$PU".fastq
