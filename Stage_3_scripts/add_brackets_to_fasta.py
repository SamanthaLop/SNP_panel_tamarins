from Bio import SeqIO
import re
import ntpath
import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
import os
import sys


# script that adds a [] around the SNP in fasta files. If running with files other than those in /proj/ \
# proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/chosen_snps_consecutive_regions/chosen_SNPs_vcfs/ 
# update the script accordingly.

# Modules to be loaded before running this script 
# module load bioinfo-tools pysam/0.16.0.1-python3.8.7 python/3.8.7 biopython/1.76-py3

#input_fasta="/proj/proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/chosen_snps_consecutive_regions/chosen_SNPs_vcfs/interval_1_chosen_snp_CM018917.1:153607702-153608404_chosen.fasta"
#output_file="/proj/proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/chosen_snps_consecutive_regions/chosen_SNPs_vcfs/interval_1_chosen_snp_CM018917.1:153607702-153608404_chosen_modified.fasta"

input_file=str(sys.argv[1])
output_file=str(sys.argv[2])


#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/species_specific_SNPs/chosen_snps_consecutive_regions/all.tsv"
#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/SNP_list_150_threshold/all.tsv"
#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/150_threshold/rm_too_close_too_end_chrom_name_org_lists/all.tsv"
#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/rerun_combine_genotype/one_vcf_per_snp/homologous_lists/all.tsv"
#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/individual_specific_SNPs/NEW_MAF/consecutive_intervals/400_threshold/all.tsv"
#snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/individual_species_specific_SNPs/lwed/lists/final_snp_list_for_brakets/all.tsv"
snp_list="/proj/proj_name/nobackup/SAM/flanking_regions/adding_to_ind_set_feb_2022/all_for_adding_brakets.tsv"

# parse snp_list, create dictionary of lists for every entry
SNPs={"chrom":[], "SNP":[], "start":[], "stop":[], "brackets_pos":[]}
with open(snp_list) as SNPlist:
    read = csv.reader(SNPlist, delimiter="\t")
    for row in read:
        SNPs["chrom"].append(row[0])
        SNPs["SNP"].append(row[1])
        SNPs["start"].append(row[2])
        SNPs["stop"].append(row[3])
        SNPs["brackets_pos"].append(int(row[1]) - int(row[2]))

# establish the start and stop positions for the input fasta
basename_input=ntpath.basename(input_file)
basename_split=re.split(r"_|-|:|\.",basename_input)
start=basename_split[6] # 6 when species-specific, 4 when individual
stop=basename_split[7] # 7 when species-scpecific, 5 when individual



for i in range(len(SNPs["SNP"])):
    start_i=SNPs["start"][i]
    stop_i=SNPs["stop"][i]
    brackets_position_i=SNPs["brackets_pos"][i]
    if start_i == start and stop_i == stop:
        fasta_sequence=SeqIO.parse(open(input_file), 'fasta')
        with open(output_file, 'w') as out_file:
            record=[]
            for fasta in fasta_sequence:
                sequence = str(fasta.seq)
                new_seqs=(str(sequence[:int(brackets_position_i)])+"["+str(sequence[brackets_position_i])+"]"+ str(sequence[int(brackets_position_i)+1:]))
                record.append(SeqRecord(Seq(str(new_seqs)), id=fasta.id))
            SeqIO.write(record, out_file, "fasta")

