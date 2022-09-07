#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J dedup_script4
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script takes merged bam file from the bwa_RG.sh script, marks duplicates and outputs a deduplicated bam file.
# the script also calls upon two other scripts, qualimap and valsam. These scripts will run qualimap and validate sam (picard)
# on the resulting deduplicated bam file. This is the final script on the automated series.

module load bioinfo-tools samtools picard

# give merged bam file as cmd input
bamfile=$1

# output in subdirectory PROCESSED_BAMS with same parent as input bam
bamdir=$(dirname $bamfile)/PROCESSED_BAMS

mkdir -p $bamdir

# grab the sampleID from the input bam
sampleID=$(basename $bamfile | cut -d "." -f 2)

# definfe reference genome 
ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna

# mark duplicates and index output
java -Xmx6g -jar "$PICARD_HOME"/picard.jar MarkDuplicates \
        I="$bamfile" \
        O="$bamdir"/dedup."$sampleID".bam \
        M="$bamdir"/metrix_"$sampleID".txt \
        ASSUME_SORT_ORDER=coordinate
samtools index "$bamdir"/dedup."$sampleID".bam

rm -f $1
rm -f mba."$SM"_"$PU".bai

# after this Axel runs two sanity checks on the processed bams. Valsam throws an error if anything is wrong with the file. Qualimap is good for
# checking depth, mapping quality etc.

sbatch -M snowy -J $(echo "$sampleID".valsam | cut -d "_" -f 1-2) --dependency=afterok:"$SLURM_JOB_ID" /proj/proj_name/nobackup/SAM/scripts/valsam.sh "$bamdir"/dedup."$sampleID".bam
sbatch -M snowy -J $(echo "$sampleID".qualimap | cut -d "_" -f 1-2) --dependency=afterok:"$SLURM_JOB_ID" /proj/proj_name/nobackup/SAM/scripts/qualimap.sh "$bamdir"/dedup."$sampleID".bam
