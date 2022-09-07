#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 5-00:00:00
#SBATCH -J qualimap
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load QualiMap

# input file is the deduplicated bam file
file=$1
name=$(basename $file)
SM=$(echo $name | cut -d _ -f 1)

# qualimap gave me errors when running, using "unset DISPLAY" fixed the issue
unset DISPLAY

# output directory to store output pdfs that come from running qualimap
OUTPUT=/proj/proj_name/nobackup/SAM/qualimap/

qualimap bamqc -nt 16 -bam $file -outdir $OUTPUT -outfile $SM.pdf -outformat pdf
