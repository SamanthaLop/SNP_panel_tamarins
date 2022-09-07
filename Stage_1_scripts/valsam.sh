#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J valsam
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

module load bioinfo-tools GATK

# use of validate sam (picard tools) reports validity of SAM or BAM files. We need out bam to be "compliant" with GATK specifications
# for downstream analysis. This tells us if the file is valid. It will not run if it isn't.

# Input file is the deduplicated bam file
file=$1
name=$(basename $file)
SM=$(echo $name | cut -d _ -f 1)

reference=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/valsam

gatk ValidateSamFile -I $file -R $reference -MODE SUMMARY -O "$out_dir"/"$SM"_valsam.txt

