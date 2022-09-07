#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -J mfeprimer
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# script to obtain a primer.fa file for mfeprimer from the batchprimer3 tsv output and run mfeprimer (indexing and in silico pcr)

output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/ANGSD_consensus_seqs/REAL_TEST/101462_round_4.2
#primer_tsv=$1
primer_filename=$(echo $(basename $1) | cut -d . -f 1)

template_fasta=$2
template_fa_filename=$(echo $(basename $2))

# using awk to get a primer list in the format that mfeprimer requires it
#awk '{print $1"\t"$12}' $primer_tsv | awk '{print ">"$1"\n"$2}' > "$output_dir"/"$primer_filename"_mfeprimer.fa

#mfeprimer_list="$output_dir"/"$primer_filename"_mfeprimer.fa
mfeprimer_list=$1

# going into the output dir:
cd $output_dir

# I had previously made a soft link, this only has to be done once:
# ln -s /home/samlope/glob/MFEprimer/mfeprimer-3.2.4-linux-amd64 mfeprimer

# index fasta to be used as template
#./mfeprimer index -i $template_fasta

# run in silico pcr
../../mfeprimer -i $mfeprimer_list -d $template_fasta > "$output_dir"/"$template_fa_filename"_"$primer_filename"_mfeoutput.txt
