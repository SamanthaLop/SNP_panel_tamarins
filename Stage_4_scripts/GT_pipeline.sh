#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 01:00:00
#SBATCH -J GT_pipeline
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# REMEMBER TO UNZIP FASTAS BEFORE USE (bgzip -d <file>) if command bgzip not found, load bioinfo-tools and then samtools

module load bioinfo-tools
module load perl/5.26.2

scripts_dir=/proj/proj_name/nobackup/SAM/scripts
input_fasta=$1
hash_dir=/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/hash
seq_test_dir=/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/seq_test
assay_txt=/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/AssayInfo_ALL.txt
probe_csv=/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/LocusInfo_ALL.csv
genotype_dir=/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/genotype

filename=$(echo $(basename $input_fasta) | cut -d . -f 1)

# HASH step
perl "$scripts_dir"/GTseq_HashSeqs.pl "$input_fasta" > "$hash_dir"/"$filename".hash
#./GTseq_HashSeqs.pl "$input_fasta" > "$hash_dir"/"$filename".hash

# seq test step
perl $scripts_dir/GTseq_SeqTest.pl $assay_txt $hash_dir/$filename.hash > $seq_test_dir/"$filename"_seqtest.csv

# genotype step
perl $scripts_dir/GTseq_Genotyper.pl $probe_csv $input_fasta > $genotype_dir/"$filename".genos

