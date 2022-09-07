#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -J merge_bams
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# Adapted from a script by Axel Jensen
module load bioinfo-tools bwa samtools GATK/4 picard

# output directory
outdir=/proj/proj_name/nobackup/SAM/tamarins_mapped

# raw reads directory
reads_dir=/proj/proj_name/nobackup/TAMARINS/sorted_tamarins

# reference genome
ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna

# scripts directory
scripts=/proj/proj_name/nobackup/SAM/scripts

# SM = Sample
SM=$(echo $(basename $1) | cut -d _ -f 1| cut -d . -f 2)

# PU = Platform unit (flowcell barcode, lane, sample barcode)
PU=$(zcat "$reads_dir"/"$SM"_1_sorted.fastq.gz | head -n 1 | cut -d ":" -f 2-4 | cut -d "@" -f 2 | sed "s|":"|"."|g")

# the next line	sends the next script to be run	AFTER this script is finished.
sbatch -J $(echo "$SM".dedup_script4 | cut -d "_" -f 1-2) --dependency=afterok:"$SLURM_JOB_ID" "$scripts"/dedup_script4.sh "$outdir"/mba."$SM"_"$PU".bam  

## script 3 4 cores
# ----- STEP 5 - merge mapped bamfile with unmapped containing readgroup data, result is mapped bamfile with all relevant metadata
java -Xmx16G -jar "$PICARD_HOME"/picard.jar MergeBamAlignment \
        R="$ref_genome" \
        UNMAPPED_BAM=$1 \
        ALIGNED_BAM=$2 \
        O="$outdir"/mba."$SM"_"$PU".bam \
        CREATE_INDEX=true \
        ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false \
        CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        ATTRIBUTES_TO_RETAIN=XS

# ----- remove unecessary intermediates
rm -f "$outdir"/mi."$SM"_"$PU".u.bam
rm -f "$outdir"/mi."$SM"_"$PU".bai
rm -f "$outdir"/sorted."$SM"_"$PU".bam
rm -f "$outdir"/sorted."$SM"_"$PU".bai
