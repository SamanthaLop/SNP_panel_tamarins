#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 30:00:00
#SBATCH -J RG_data_script1
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# adapted from a script provided by Axel Jensen
# This script -if needed- indexes and/or creates the dictionary of the reference genome. It converts the fastq files into
# sam files to which the read group data is added. It also marks the ilumina adapters and reverts the sam file back to fastq format
# to be able to map with bwa. The 1 at the end of the script name indicates this is the first script in a series of 4
# automated scripts.

module load bioinfo-tools GATK/4 picard

# script directory
scripts=/proj/proj_name/nobackup/SAM/scripts

# the two input files are the paired end fastq files for each sample
raw_read_1=$1
raw_read_2=$2     # from axels script "raw_read_2=${1/_1_sorted.fastq.gz/_2_sorted.fastw.gz}

# output directory
outdir=/proj/proj_name/nobackup/SAM/tamarins_mapped

# reference genome
ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna

# index reference genome, in this case this is already done so no need to run it, if you need to run this, just delete the "#"
#bwa index $ref_genome

# create dictionary of reference genome. This just has to be done once.
#java -jar "$PICARD_HOME"/picard.jar CreateSequenceDictionary \
#       R=$ref_genome \
#       O=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.dict

# read group data
name=$(basename $raw_read_1)

# SM = Sample
SM=$(echo $name | cut -d _ -f 1)
# first file, gives me this: 101455

# ID = Read group identifier
ID=$(zcat "$raw_read_1" | head -n 1 | cut -d ":" -f 1,2-4 | cut -d "@" -f 2 | sed "s|":"|"."|g")
# first file, gives me this: A00893.127.HTKKYDSXY.4

# PU = Platform unit (flowcell barcode, lane, sample barcode)
PU=$(zcat "$raw_read_1" | head -n 1 | cut -d ":" -f 2-4 | cut -d "@" -f 2 | sed "s|":"|"."|g")
# first file, gives me this: 127.HTKKYDSXY.4

# PL = Platform/technology used to produce the read
PL=ILLUMINA

# LB = DNA preparation library identifier
LB="$SM"

# PI = insert size, Axel: I think this value should work for your data as well.
PI=330

# the next line sends the next script to be run AFTER this script is finished.
sbatch -J $(echo "$SM".bwa_script2) --dependency=afterok:"$SLURM_JOB_ID" "$scripts"/bwa_script2.sh "$outdir"/"$SM"_"$PU".fastq

# ----- STEP 1 - convert the fastw files to a sam and provide read group information
java -Xmx6G -Xms6G -jar "$PICARD_HOME"/picard.jar FastqToSam \
        F1=$raw_read_1 \
        F2=$raw_read_2 \
        SM="$SM" \
        RG="$ID" \
        PU="$PU" \
        PL="$PL" \
        LB="$LB" \
        PI="$PI" \
        O=$outdir/"$SM"_"$PU".bam

# ----- STEP 2 - mark illumina adapter content on this samfile, output is ubam (unmapped bam)
java -Xmx6G -Xms6G -jar "$PICARD_HOME"/picard.jar MarkIlluminaAdapters \
        I=$outdir/"$SM"_"$PU".bam \
        O="$outdir"/mi."$SM"_"$PU".u.bam \
        M="$outdir"/mi."$SM"_"$PU".metrix.txt

# ----- remove unnecessary intermediate
rm -r -f $outdir/"$SM"_"$PU".bam

# ----- STEP 3 - revert sam file back to interleaved fastq to map in bwa
java -Xmx6G -Xms6G -jar "$PICARD_HOME"/picard.jar SamToFastq \
        I="$outdir"/mi."$SM"_"$PU".u.bam \
        FASTQ="$outdir"/"$SM"_"$PU".fastq \
        CLIPPING_ATTRIBUTE=XT \
        CLIPPING_ACTION=2 \
        INTERLEAVE=true \
        NON_PF=true
