# hethom, homhom
import csv
from email import header
import glob
from os import spawnlpe
from pickletools import int4
from traceback import format_exception_only
import pandas as pd
from pandas import DataFrame
from collections import defaultdict
from Bio import SeqIO
import re

dict_of_inds = {}

for file in glob.glob("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/updated_geno/*.genos"):
    df = pd.read_csv(file, sep = "\t")
    ind_info = '_'.join(file.split('/')[-1].split('_')[0:2])
    dict_of_snps = {}
    for index, row in df.iterrows():
        chrom = row["Chrom"]
        pos = row["Position"]
        bp3_index = row["bp3_index"]
        a1 = row["A1counts"].split('=')[0]
        a1counts = row["A1counts"].split('=')[1]
        a2 = row["A2counts"].split('=')[0]
        a2counts = row["A2counts"].split('=')[1]
        genotype = row["Genotype"]
        genotype_class = row["Genotype_class"]
        if bp3_index not in dict_of_snps:
            dict_of_snps[bp3_index] = a1, a1counts, a2, a2counts, genotype, genotype_class
    if ind_info not in dict_of_inds:
        dict_of_inds[ind_info] = dict_of_snps


ind_id_dual = [21,22,23,24,25,26,29,30,31,32,34,37,38,39,42,43,44,46,47,49,51,53,54,55,58,60,61,62,64,65,67,69,70,74,75,76,78,80,81,82,83,86,88,89,90,91,92,94,96,98,100,101,103,104,105,107,108,109,110,111,112,113,114,117,118,119,120,123,124,128,129,130,131,132,133,134,135,137,138,139,140,141,142,144,145,146,147,148,149,150,152,153,154,156,157,158,159,160,161,162,164,165,166,167,168,169,171,173,174,176,177,178,179,181,183,184,187,188,189,190,192,193,358,360,361,363,364,367,369,372,373,376,378,379,380,381,382,383,384,385,388,389,390,392,393,394,395,396,399,400,401]
ind_id_lwed = [226,227,228,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,248,249,250,251,252,253,254,255,256,257,258,259,260,261,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291]
ind_id_simp = [292,293,295,296,297,299,300,302,303,304,306,307,308,310,311,312,313,315,316,317,318,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,341,343,344,345,346,347,348,349,352,353,355,356,357]

failed_ind_id = []
failed_lwed = []
failed_simp = []

for snp in no_reads_snps:
    if snp in ind_id_dual:
        failed_ind_id.append(snp)
    if snp in ind_id_lwed:
        failed_lwed.append(snp)
    if snp in ind_id_simp:
        failed_simp.append(snp)

primers_file="/proj/proj_name/nobackup/SAM/flanking_regions/final_merged_fastas/19.02.2022_fastas_for_merging/simulating_seqs_to_troubleshoot_gt_seq/primers.tsv"
fasta_new_snps="/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/concatenated_ALL_19.02.2022_adding_brakets_replace_snps_updated4.6.22.fasta"

primers = pd.read_csv(primers_file, sep = '\t')

def reverse_complement(dna):
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join([complement[base] for base in dna[::-1]])

correct = 0
incorrect = 0

snps_that_arent_working =[]

with open(fasta_new_snps) as f:
    SeqDict = {}
    for record in SeqIO.parse(f, "fasta"):
        sequence = str(record.seq)
        snp_name = str(record.id)
        bp3_index = snp_name.split('_')[-1]
        for index, row in primers.iterrows():
            forward_p=str(row["forward"])
            reverse_p=str(reverse_complement(row["reverse"]))
            if forward_p in sequence and reverse_p in sequence:
                amplicon = re.search('{}(.*?){}'.format(forward_p, reverse_p), sequence).group()
                if '[' in amplicon and ']' in amplicon:
                    correct = correct + 1
                else:
                    incorrect = incorrect + 1
                    snps_that_arent_working.append(bp3_index)


bp3_err_ind_id = []
bp3_err_lwed = []
bp3_err_simp = []

for snp in snps_that_arent_working:
    if int(snp) in ind_id_dual:
        bp3_err_ind_id.append(snp)
    if int(snp) in ind_id_lwed:
        bp3_err_lwed.append(snp)
    if int(snp) in ind_id_simp:
        bp3_err_simp.append(snp)


under_performin_snps = [338, 361, 138, 67, 42, 37, 26, 280, 326, 174]

dual_list = [] # length = 145   NOW 146  # NOW 139
for snp in ind_id_dual:
    if str(snp) not in snps_that_arent_working and snp not in failed_ind_id and snp not in under_performin_snps:
        dual_list.append(snp)

lwed_list = [] # length = 25   NOW 37 # now 37
for snp in ind_id_lwed:
    if str(snp) not in snps_that_arent_working and snp not in failed_lwed and snp not in under_performin_snps:
        lwed_list.append(snp)

simp_list = [] # length = 19    # NOW 41 # now 40
for snp in ind_id_simp:
    if str(snp) not in snps_that_arent_working and snp not in failed_simp and snp not in under_performin_snps:
        simp_list.append(snp)

ind_spp= {'tamOpt-001_S1': "SIMP", "tamOpt-016_S14": "SIMP", "tamOpt-008_S8": "LWED", "tamOpt-018_S17": "SIMP", "tamOpt-025_S16": "NA",\
            "tamOpt-032_S30": "LWED", "tamOpt-027_S20": "LWED", "tamOpt-026_S18": "LWED", "tamOpt-017_S15": "SIMP","tamOpt-031_S28": "LWED", \
                "tamOpt-002_S2": "LWED", "tamOpt-020_S21": "SIMP", "tamOpt-004_S4": "NA", "tamOpt-029_S24": "SIMP", "tamOpt-010_S9": "LWED",\
                    "tamOpt-023_S27": "SIMP", "tamOpt-015_S13": "SIMP", "tamOpt-003_S3": "LWED", "tamOpt-033_S31": "LWED", "tamOpt-006_S6": "SIMP", \
                        "tamOpt-012_S10": "LWED", "tamOpt-019_S19": "SIMP", "tamOpt-034_S32": "SIMP", "tamOpt-030_S26": "LWED", "tamOpt-024_S29": "SIMP",\
                            "tamOpt-022_S25": "NA", "tamOpt-007_S7": "NA", "tamOpt-021_S23": "SIMP", "tamOpt-013_S11": "LWED", "tamOpt-005_S5": "LWED",\
                                "tamOpt-014_S12": "LWED","tamOpt-028_S22": "SIMP","tamOpt-035_S33": "LWED"}

ind_gen = {}
snp_pos_lists = {}


for ind in dict_of_inds:
    ind_list = []
    snp_pos = []
    spp = ind_spp[ind]
    for snp in dict_of_inds[ind]:
        if (snp in dual_list) or (snp in lwed_list) or (snp in simp_list):
            #print(dict_of_inds[ind][snp])
            if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 20:
                a1_allele = dict_of_inds[ind][snp][0]
                a2_allele = dict_of_inds[ind][snp][2]
                genotype_class = dict_of_inds[ind][snp][5]
                if genotype_class == "A1HOM":
                    ind_list.append(a1_allele)
                    snp_pos.append(snp)
                elif genotype_class == "A2HOM":
                    ind_list.append(a2_allele)
                    snp_pos.append(snp)
                elif genotype_class == "HET":
                    if (a1_allele == "C" and a2_allele == "A") or (a1_allele == "A" and a2_allele == "C"):
                        ind_list.append("M")
                        snp_pos.append(snp)
                    elif (a1_allele == "A" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "A"):
                        ind_list.append("R")
                        snp_pos.append(snp)
                    elif (a1_allele == "A" and a2_allele == "T") or (a1_allele == "T" and a2_allele == "A"):
                        ind_list.append("W")
                        snp_pos.append(snp)
                    elif (a1_allele == "C" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "C"):
                        ind_list.append("S")
                        snp_pos.append(snp)
                    elif (a1_allele == "C" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "C"):
                        ind_list.append("Y")
                        snp_pos.append(snp)
                    elif (a1_allele == "G" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "G"):
                        ind_list.append("K")
                        snp_pos.append(snp)
                    else:
                        ind_list.append("N")
                        snp_pos.append(snp) # multivariate site
                else:
                    ind_list.append("N") # nan genotype change to x recalculate matrix
                    snp_pos.append(snp)
            elif (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) < 20:
                ind_list.append("N") # less than 10 reads
                snp_pos.append(snp)
    snp_pos_lists[ind] = snp_pos
    ind_gen[ind]=ind_list

ind_gen_joined = {}
for ind in ind_gen:
    string = ''.join(ind_gen[ind])
    ind_gen_joined[ind] = string

output = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/dual_ind_test_1.tsv"

with open(output, "w") as out:
    out.write("sample" + "\t" + "dual_snps_genotype" + "\n")
    for ind in ind_gen_joined:
        out.write(ind + "\t" + ind_gen_joined[ind] + "\n")

def compare(chain1, chain2):
    dist_hom_vs_het = 0
    dist_hom_vs_hom = 0
    dist_extra = 0
    for i in range(len(chain1)):
        if (chain1[i] == "X") or (chain1[i] == "N") or (chain2[i] == "X") or (chain2[i] == "N"):
            dist_hom_vs_het = dist_hom_vs_het + 0
            dist_hom_vs_hom = dist_hom_vs_hom + 0
        else:
            if chain1[i] in ["A", "C", "T", "G"]:
                if chain2[i] in ["A", "C", "T", "G"]:
                    if chain1[i] == chain2[i]:
                        dist_hom_vs_hom = dist_hom_vs_hom + 0
                    else:
                        dist_hom_vs_hom = dist_hom_vs_hom + 2
                elif chain2[i] in ["R", "Y", "S", "W", "K", "M"]:
                    if chain2[i] == "R":
                        if chain1[i] in ["A", "G"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["C", "T"]:
                            dist_extra = dist_extra + 2
                    if chain2[i] == "Y":
                        if chain1[i] in ["C", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["A", "G"]:
                            dist_extra = dist_extra + 2
                    if chain2[i] == "S":
                        if chain1[i] in ["G", "C"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["T", "A"]:
                            dist_extra = dist_extra + 2
                    if chain2[i] == "W":
                        if chain1[i] in ["A", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["C", "G"]:
                            dist_extra = dist_extra + 2
                    if chain2[i] == "K":
                        if chain1[i] in ["G", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["C", "A"]:
                            dist_extra = dist_extra + 2
                    if chain2[i] == "M":
                        if chain1[i] in ["A", "C"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain1[i] in ["T", "G"]:
                            dist_extra = dist_extra + 2
            elif chain1[i] in ["R", "Y", "S", "W", "K", "M"]:
                if chain2[i] in ["R", "Y", "S", "W", "K", "M"]:
                    if chain1[i] == chain2[i]:
                        dist_extra = dist_extra + 0
                    else:
                        dist_extra = dist_extra + 2
                elif chain2[i] in ["A", "C", "T", "G"]:
                    if chain1[i] == "R":
                        if chain2[i] in ["A", "G"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain2[i] in ["T", "C"]:
                            dist_extra = dist_extra + 2
                    if chain1[i] == "Y":
                        if chain2[i] in ["C", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain2[i] in ["A", "G"]:
                            dist_extra = dist_extra + 2
                    if chain1[i] == "S":
                        if chain2[i] in ["G", "C"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain2[i] in ["A", "T"]:
                            dist_extra = dist_extra + 2
                    if chain1[i] == "W":
                        if chain2[i] in ["A", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain2[i] in ["C", "G"]:
                            dist_extra = dist_extra + 2
                    if chain1[i] == "K":
                        if chain2[i] in ["G", "T"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        elif chain2[i] in ["C", "A"]:
                            dist_extra = dist_extra + 2
                    if chain1[i] == "M":
                        if chain2[i] in ["A", "C"]:
                            dist_hom_vs_het = dist_hom_vs_het + 1
                        if chain2[i] in ["G", "T"]:
                            dist_extra = dist_extra + 2
    return dist_extra, dist_hom_vs_het, dist_hom_vs_hom

############## LWED
ind_gen_lwed = {}
snp_pos_lists = {}

for ind in dict_of_inds:
    ind_list = []
    snp_pos = []
    spp = ind_spp[ind]
    if ind_spp[ind] == "LWED":
        for snp in dict_of_inds[ind]:
            if (snp in dual_list) or (snp in lwed_list): #or (snp in simp_list):
                #print(dict_of_inds[ind][snp])
                if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 10:
                    a1_allele = dict_of_inds[ind][snp][0]
                    a2_allele = dict_of_inds[ind][snp][2]
                    genotype_class = dict_of_inds[ind][snp][5]
                    if genotype_class == "A1HOM":
                        ind_list.append(a1_allele)
                        snp_pos.append(snp)
                    elif genotype_class == "A2HOM":
                        ind_list.append(a2_allele)
                        snp_pos.append(snp)
                    elif genotype_class == "HET":
                        if (a1_allele == "C" and a2_allele == "A") or (a1_allele == "A" and a2_allele == "C"):
                            ind_list.append("M")
                            snp_pos.append(snp)
                        elif (a1_allele == "A" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "A"):
                            ind_list.append("R")
                            snp_pos.append(snp)
                        elif (a1_allele == "A" and a2_allele == "T") or (a1_allele == "T" and a2_allele == "A"):
                            ind_list.append("W")
                            snp_pos.append(snp)
                        elif (a1_allele == "C" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "C"):
                            ind_list.append("S")
                            snp_pos.append(snp)
                        elif (a1_allele == "C" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "C"):
                            ind_list.append("Y")
                            snp_pos.append(snp)
                        elif (a1_allele == "G" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "G"):
                            ind_list.append("K")
                            snp_pos.append(snp)
                        else:
                            ind_list.append("N")
                            snp_pos.append(snp) # multivariate site
                    else:
                        ind_list.append("N") # nan genotype change to x recalculate matrix
                        snp_pos.append(snp)
                elif (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) < 10:
                    ind_list.append("N") # less than 10 reads
                    snp_pos.append(snp)
        snp_pos_lists[ind] = snp_pos
        ind_gen_lwed[ind]=ind_list

# LWED and dual snps result in sequences of 176 in length:
#>>> snp_pos
#[21, 22, 24, 25, 29, 30, 31, 32, 34, 38, 39, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 69, 70, 74, 75, 76, 78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 176, 177, 178, 179, 181, 183, 184, 187, 188, 189, 190, 192, 193, 358, 360, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394, 395, 396, 399, 400, 401, 226, 232, 233, 235, 237, 239, 240, 244, 246, 248, 249, 250, 254, 255, 257, 259, 263, 264, 266, 268, 270, 271, 272, 274, 275, 276, 278, 279, 281, 282, 283, 284, 285, 286, 288, 289, 290]

under_performing_inds = ["tamOpt-035_S33", "tamOpt-007_S7", 'tamOpt-022_S25', 'tamOpt-025_S16', "tamOpt-004_S4"]

ind_gen_joined_lwed = {}
for ind in ind_gen_lwed:
    if ind not in under_performing_inds:
        string = ''.join(ind_gen_lwed[ind])
        ind_gen_joined_lwed[ind] = string

output = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/custom_snpsseqs_LWED.tsv"
with open(output, "w") as out:
    out.write("sample" + "\t" + "dual_snps_genotype" + "\n")
    for ind in ind_gen_joined_lwed:
        out.write(ind + "\t" + ind_gen_joined_lwed[ind] + "\n")

individuals_lwed = []
for i in ind_gen_joined_lwed:
    individuals_lwed.append(i)

# as a dictionary
ind_genotypes_dual_dict_lwed_extra = {}
ind_genotypes_dual_dict_lwed_hethom = {}
ind_genotypes_dual_dict_lwed_homhom = {}
for i in ind_gen_joined_lwed:
    ind_gens_extra = []
    ind_gens_hethom = []
    ind_gens_homhom = []
    seq1 = ind_gen_joined_lwed[i]
    for j in ind_gen_joined_lwed:
        seq2 = ind_gen_joined_lwed[j]
        d_extra = compare(seq1, seq2)[0]
        d_hethom = compare(seq1,seq2)[1]
        d_homhom = compare(seq1,seq2)[2]
        ind_gens_extra.append(str(d_extra))
        ind_gens_hethom.append(str(d_hethom))
        ind_gens_homhom.append(str(d_homhom))
    ind_genotypes_dual_dict_lwed_extra[i] = ind_gens_extra
    ind_genotypes_dual_dict_lwed_hethom[i] = ind_gens_hethom
    ind_genotypes_dual_dict_lwed_homhom[i] = ind_gens_homhom

# extra is empty so let's make two matrices now.

geno_dual_lwed_hethom = pd.DataFrame.from_dict(ind_genotypes_dual_dict_lwed_hethom)
geno_dual_lwed_homhom = pd.DataFrame.from_dict(ind_genotypes_dual_dict_lwed_homhom)

t_lwed_hethom = geno_dual_lwed_hethom.transpose().astype(int)
t_lwed_homhom = geno_dual_lwed_homhom.transpose().astype(int)

t_lwed_hethom.columns = individuals_lwed
t_lwed_homhom.columns = individuals_lwed


new_t_lwed_hethom = t_lwed_hethom.rename(columns = {'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"}, index={'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"})
new_t_lwed_homhom = t_lwed_homhom.rename(columns = {'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"}, index={'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"})

import numpy as np
lwed_mask_hethom = np.zeros_like(new_t_lwed_hethom, dtype=np.bool)
lwed_mask_homhom = np.zeros_like(new_t_lwed_homhom, dtype=np.bool)


lwed_mask_hethom[np.triu_indices_from(lwed_mask_hethom)] = True
lwed_mask_homhom[np.triu_indices_from(lwed_mask_homhom)] = True

lwed_mask_hethom[np.diag_indices_from(lwed_mask_hethom)] = False
lwed_mask_homhom[np.diag_indices_from(lwed_mask_homhom)] = False

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_lwed_hethom, mask=lwed_mask_hethom, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for L. weddelli - Dual and lwed SNP sets - hethom")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_LWED_duallwedsets_hethom_20threshold.png", dpi=300)


sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_lwed_homhom, mask=lwed_mask_homhom, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for L. weddelli - Dual and lwed SNP sets - homhom")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_LWED_duallwedsets_homhom_20threshold.png", dpi=300)

############## SIMP and dual SNP sets:
ind_gen_simp = {}
snp_pos_lists = {}

for ind in dict_of_inds:
    ind_list = []
    snp_pos = []
    spp = ind_spp[ind]
    if ind_spp[ind] == "SIMP":
        for snp in dict_of_inds[ind]:
            if (snp in dual_list) or (snp in simp_list): #or (snp in simp_list):
                #print(dict_of_inds[ind][snp])
                if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 10:
                    a1_allele = dict_of_inds[ind][snp][0]
                    a2_allele = dict_of_inds[ind][snp][2]
                    genotype_class = dict_of_inds[ind][snp][5]
                    if genotype_class == "A1HOM":
                        ind_list.append(a1_allele)
                        snp_pos.append(snp)
                    elif genotype_class == "A2HOM":
                        ind_list.append(a2_allele)
                        snp_pos.append(snp)
                    elif genotype_class == "HET":
                        if (a1_allele == "C" and a2_allele == "A") or (a1_allele == "A" and a2_allele == "C"):
                            ind_list.append("M")
                            snp_pos.append(snp)
                        elif (a1_allele == "A" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "A"):
                            ind_list.append("R")
                            snp_pos.append(snp)
                        elif (a1_allele == "A" and a2_allele == "T") or (a1_allele == "T" and a2_allele == "A"):
                            ind_list.append("W")
                            snp_pos.append(snp)
                        elif (a1_allele == "C" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "C"):
                            ind_list.append("S")
                            snp_pos.append(snp)
                        elif (a1_allele == "C" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "C"):
                            ind_list.append("Y")
                            snp_pos.append(snp)
                        elif (a1_allele == "G" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "G"):
                            ind_list.append("K")
                            snp_pos.append(snp)
                        else:
                            ind_list.append("N")
                            snp_pos.append(snp) # multivariate site
                    else:
                        ind_list.append("N") # nan genotype change to x recalculate matrix
                        snp_pos.append(snp)
                elif (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) < 10:
                    ind_list.append("N") # less than 10 reads
                    snp_pos.append(snp)
        snp_pos_lists[ind] = snp_pos
        ind_gen_simp[ind]=ind_list

# snp_pos => [21, 22, 24, 25, 29, 30, 31, 32, 34, 38, 39, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 69, 70, 74, 75, 76, 78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 176, 177, 178, 179, 181, 183, 184, 187, 188, 189, 190, 192, 193, 358, 360, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394, 395, 396, 399, 400, 401, 293, 296, 297, 299, 302, 304, 306, 307, 308, 311, 312, 313, 315, 317, 320, 321, 322, 323, 324, 325, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 341, 345, 346, 348, 352, 355, 356, 357]

ind_gen_joined_simp = {}
for ind in ind_gen_simp:
    if ind not in under_performing_inds:
        string = ''.join(ind_gen_simp[ind])
        ind_gen_joined_simp[ind] = string

output = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/custom_snpsseqs_SIMP.tsv"
with open(output, "w") as out:
    out.write("sample" + "\t" + "dual_snps_genotype" + "\n")
    for ind in ind_gen_joined_simp:
        out.write(ind + "\t" + ind_gen_joined_simp[ind] + "\n")

individuals_simp = []
for i in ind_gen_joined_simp:
    individuals_simp.append(i)


# as a dictionary
ind_genotypes_dual_dict_simp_extra = {}
ind_genotypes_dual_dict_simp_hethom = {}
ind_genotypes_dual_dict_simp_homhom = {}
for i in ind_gen_joined_simp:
    ind_gens_extra = []
    ind_gens_hethom = []
    ind_gens_homhom = []
    seq1 = ind_gen_joined_simp[i]
    for j in ind_gen_joined_simp:
        seq2 = ind_gen_joined_simp[j]
        d_extra = compare(seq1, seq2)[0]
        d_hethom = compare(seq1,seq2)[1]
        d_homhom = compare(seq1,seq2)[2]
        ind_gens_extra.append(str(d_extra))
        ind_gens_hethom.append(str(d_hethom))
        ind_gens_homhom.append(str(d_homhom))
    ind_genotypes_dual_dict_simp_extra[i] = ind_gens_extra
    ind_genotypes_dual_dict_simp_hethom[i] = ind_gens_hethom
    ind_genotypes_dual_dict_simp_homhom[i] = ind_gens_homhom

# extra is empty so let's make two matrices now.

geno_dual_simp_hethom = pd.DataFrame.from_dict(ind_genotypes_dual_dict_simp_hethom)
geno_dual_simp_homhom = pd.DataFrame.from_dict(ind_genotypes_dual_dict_simp_homhom)

t_simp_hethom = geno_dual_simp_hethom.transpose().astype(int)
t_simp_homhom = geno_dual_simp_homhom.transpose().astype(int)

t_simp_hethom.columns = individuals_simp
t_simp_homhom.columns = individuals_simp

new_t_simp_hethom = t_simp_hethom.rename(columns = {'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"}, index={'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"})
new_t_simp_homhom = t_simp_homhom.rename(columns = {'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"}, index={'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"})


import numpy as np
simp_mask_hethom = np.zeros_like(new_t_simp_hethom, dtype=np.bool)
simp_mask_homhom = np.zeros_like(new_t_simp_homhom, dtype=np.bool)

simp_mask_hethom[np.triu_indices_from(simp_mask_hethom)] = True
simp_mask_homhom[np.triu_indices_from(simp_mask_homhom)] = True

simp_mask_hethom[np.diag_indices_from(simp_mask_hethom)] = False
simp_mask_homhom[np.diag_indices_from(simp_mask_homhom)] = False

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_simp_hethom, mask=simp_mask_hethom, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for S. imperator - Dual and simp SNP sets - hethom")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_SIMP_duallwedsets_hethom_20threshold.png", dpi=300)


sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_simp_homhom, mask=simp_mask_homhom, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for S. imperator - Dual and simp SNP sets - homhom")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_SIMP_duallwedsets_homhom_20threshold.png", dpi=300)

# plotting both matrices in same one
simp_mask_hethom = np.zeros_like(new_t_simp_hethom, dtype=np.bool)
simp_mask_homhom = np.zeros_like(new_t_simp_homhom, dtype=np.bool)

simp_mask_hethom[np.triu_indices_from(simp_mask_hethom)] = True
simp_mask_homhom[np.tril_indices_from(simp_mask_homhom)] = True

simp_mask_hethom[np.diag_indices_from(simp_mask_hethom)] = False
simp_mask_homhom[np.diag_indices_from(simp_mask_homhom)] = False

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
#sns.set(rc = {"figure.figsize":(10,8)})
fig, ax = plt.subplots(figsize=(11,14))
ax = sns.heatmap(new_t_simp_hethom, mask=simp_mask_hethom, annot= True, annot_kws={"size":11}, cmap="Blues", square = True, cbar_kws={"label": "Heterozygous vs. homozygous site distance points (Possible allelic dropout)", "orientation": "horizontal", "shrink": 0.4, "pad": 0.009})
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
ax.tick_params(axis="y", labelrotation = 0)
ax = sns.heatmap(new_t_simp_homhom, mask=simp_mask_homhom, annot= True, annot_kws={"size":11}, cmap="Greys", square = True, cbar_kws={"label": "Homozygous vs. homozygous site distance points (Possible true allelic difference)", "orientation": "horizontal", "shrink": 0.4, "pad": 0.09})
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
ax.tick_params(axis="y", labelrotation = 0)
ax.set_title("Distance matrices of concatenated indSNP sets for S. imperator \n Dual and SIMP indSNPs")
#ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_SIMP_duallwedsets_hethom_homhom.png", dpi=300)

#LWED
lwed_mask_hethom = np.zeros_like(new_t_lwed_hethom, dtype=np.bool)
lwed_mask_homhom = np.zeros_like(new_t_lwed_homhom, dtype=np.bool)

lwed_mask_hethom[np.triu_indices_from(lwed_mask_hethom)] = True
lwed_mask_homhom[np.tril_indices_from(lwed_mask_homhom)] = True

lwed_mask_hethom[np.diag_indices_from(lwed_mask_hethom)] = False
lwed_mask_homhom[np.diag_indices_from(lwed_mask_homhom)] = False

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
#sns.set(rc = {"figure.figsize":(10,8)})
fig, ax = plt.subplots(figsize=(11,14))
ax = sns.heatmap(new_t_lwed_hethom, mask=lwed_mask_hethom, annot= True, annot_kws={"size":11}, cmap="Blues", square = True, cbar_kws={"label": "Heterozygous vs. homozygous site distance points (Possible allelic dropout)", "orientation": "horizontal", "shrink": 0.4, "pad": 0.009})
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
ax.tick_params(axis="y", labelrotation = 0)
ax = sns.heatmap(new_t_lwed_homhom, mask=lwed_mask_homhom, annot= True, annot_kws={"size":11}, cmap="Greys", square = True, cbar_kws={"label": "Homozygous vs. homozygous site distance points (Possible true allelic difference)", "orientation": "horizontal", "shrink": 0.4, "pad": 0.09})
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
ax.tick_params(axis="y", labelrotation = 0)
ax.set_title("Distance matrices of concatenated indSNP sets for L. weddelli \n Dual and LWED indSNPs")
#ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_LWED_duallwedsets_hethom_homhom.png", dpi=300)

