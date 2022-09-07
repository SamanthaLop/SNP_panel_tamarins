# ind id snps
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

print(glob.glob("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/genotype/updated_geno/*.genos"))

# make the dictionary with the genotype information

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

# transpose it to be snp and then ind

snp_dir = DataFrame(dict_of_inds).transpose().to_dict()

# let's select the snps that are producing no reads

count_dict = {}
no_reads_snps = []

for snp in snp_dir:
    failed_count_zero = 0
    no_reads_list = []
    for ind in snp_dir[snp]:
        if int(snp_dir[snp][ind][1]) + int(snp_dir[snp][ind][3]) == 0:
            failed_count_zero = failed_count_zero + 1
            no_reads_list.append(snp)
    if len(no_reads_list) == 33:
        no_reads_snps.append(snp)


# this is the list of snps not producing reads:
#no_reads_snps = [3, 6, 7, 10, 23, 51, 58, 96, 287]

# let's see which of these are ind id snps:

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

# ind id we have 4 SNPs:
#>>> failed_ind_id
#[23, 51, 58, 96]

# lwed ind id we have 1 SNP:
#>>> failed_lwed
#[287]

# simp ind id we have 0 SNPs:
#>>> failed_simp
#[]

# we also need the list of snps that had the bp3 error:

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

# we have the list now of bp3 failing snps:
#>>> snps_that_arent_working
#['390', '393', '201', '216', '220', '292', '293', '295', '296', '297', '300', '303', '307', '308', '310', \
#    '311', '312', '315', '316', '317', '318', '322', '323', '325', '329', '330', '331', '334', '335', '338', \
#        '341', '343', '344', '345', '347', '349', '352', '353', '356', '226', '227', '228', '230', '231', '234', \
#            '236', '240', '241', '242', '243', '244', '245', '246', '248', '250', '251', '252', '253', '255', '256', \
#                '258', '260', '261', '263', '265', '266', '267', '269', '272', '273', '275', '277', '279', '280', '283', '291']

# NOW
# >>> snps_that_arent_working
['390', '201', '216', '220', '292', '295', '300', '303', '310', '330', '338', '343', '344', '347', '349', \
    '353', '227', '228', '230', '231', '234', '236', '241', '242', '243', '245', '251', '252', '253', '256', \
        '258', '260', '261', '265', '267', '269', '273', '277', '280', '291']





# let's get three lists for ind id, with snps that are both producing reads, and that aren't in the not working list:


dual_list = [] # length = 145   NOW 146
for snp in ind_id_dual:
    if str(snp) not in snps_that_arent_working and snp not in failed_ind_id:
        dual_list.append(snp)

lwed_list = [] # length = 25   NOW 37
for snp in ind_id_lwed:
    if str(snp) not in snps_that_arent_working and snp not in failed_lwed:
        lwed_list.append(snp)

simp_list = [] # length = 19    # NOW 41
for snp in ind_id_simp:
    if str(snp) not in snps_that_arent_working and snp not in failed_simp:
        simp_list.append(snp)

# But we also have some variants that are producing less reads than we want (<5 reads in <5 individuals):
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

# these are the snps we are going to use for ind id. 
# we now know the inferred species:
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
            if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 20: # changed from 10 to 20
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
    ind_gen[ind]=ind_list

# wee test to make an lwed heatmap and a simp one:
ind_gen = {}

for ind in dict_of_inds:
    ind_list = []
    spp = ind_spp[ind]
    if spp == "LWED" or spp == "NA": # or SIMP
        for snp in dict_of_inds[ind]:
            if snp in dual_list or snp in lwed_list: #or snp in simp_list:
                #print(dict_of_inds[ind][snp])
                if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 10:
                    a1_allele = dict_of_inds[ind][snp][0]
                    a2_allele = dict_of_inds[ind][snp][2]
                    genotype_class = dict_of_inds[ind][snp][5]
                    if genotype_class == "A1HOM":
                        ind_list.append(a1_allele)
                    elif genotype_class == "A2HOM":
                        ind_list.append(a2_allele)
                    elif genotype_class == "HET":
                        if (a1_allele == "C" and a2_allele == "A") or (a1_allele == "A" and a2_allele == "C"):
                            ind_list.append("M")
                        elif (a1_allele == "A" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "A"):
                            ind_list.append("R")
                        elif (a1_allele == "A" and a2_allele == "T") or (a1_allele == "T" and a2_allele == "A"):
                            ind_list.append("W")
                        elif (a1_allele == "C" and a2_allele == "G") or (a1_allele == "G" and a2_allele == "C"):
                            ind_list.append("S")
                        elif (a1_allele == "C" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "C"):
                            ind_list.append("Y")
                        elif (a1_allele == "G" and a2_allele == "T") or (a2_allele == "T" and a2_allele == "G"):
                            ind_list.append("K")
                        else:
                            ind_list.append("N") # multivariate site
                    else:
                        ind_list.append("X") # nan genotype change to x recalculate matrix
                elif (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) < 10:
                    ind_list.append("X") # less than 10 reads
        ind_gen[ind]=ind_list

# the sequences contain these SNPs in this order:
# [21, 22, 24, 25, 29, 30, 31, 32, 34, 38, 39, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 69, 70, 74, 75, 76, 78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 176, 177, 178, 179, 181, 183, 184, 187, 188, 189, 190, 192, 193, 358, 360, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394, 395, 396, 399, 400, 401, 293, 296, 297, 299, 302, 304, 306, 307, 308, 311, 312, 313, 315, 317, 320, 321, 322, 323, 324, 325, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 341, 345, 346, 348, 352, 355, 356, 357, 226, 232, 233, 235, 237, 239, 240, 244, 246, 248, 249, 250, 254, 255, 257, 259, 263, 264, 266, 268, 270, 271, 272, 274, 275, 276, 278, 279, 281, 282, 283, 284, 285, 286, 288, 289, 290]
ind_gen_joined = {}
for ind in ind_gen:
    string = ''.join(ind_gen[ind])
    ind_gen_joined[ind] = string

output = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/dual_ind_test_1_20threshold.tsv"

with open(output, "w") as out:
    out.write("sample" + "\t" + "dual_snps_genotype" + "\n")
    for ind in ind_gen_joined:
        out.write(ind + "\t" + ind_gen_joined[ind] + "\n")

# okay let's make a function that compares two strings:
def compare(chain1, chain2):
    dist = 0
    for i in range(len(chain1)):
        if (chain1[i] == "X") or (chain1[i] == "N") or (chain2[i] == "X") or (chain2[i] == "N"):
            dist = dist + 0
        else:
            if chain1[i] == chain2[i]:
                dist = dist + 0
            elif chain1[i] != chain2[i]:
                dist = dist + 1
    return dist

# let's make a better function that accounts for allelic drop out!  _________________NEW FUNCTION___________
def compare(chain1, chain2):
    dist = 0
    for i in range(len(chain1)):
        if (chain1[i] == "X") or (chain1[i] == "N") or (chain2[i] == "X") or (chain2[i] == "N"):
            dist = dist + 0
        else:
            if chain1[i] in ["A", "C", "T", "G"]:
                if chain2[i] in ["A", "C", "T", "G"]:
                    if chain1[i] == chain2[i]:
                        dist = dist + 0
                    else:
                        dist = dist + 2
                elif chain2[i] in ["R", "Y", "S", "W", "K", "M"]:
                    if chain2[i] == "R":
                        if chain1[i] in ["A", "G"]:
                            dist = dist + 1
                        elif chain1[i] in ["C", "T"]:
                            dist = dist + 2
                    if chain2[i] == "Y":
                        if chain1[i] in ["C", "T"]:
                            dist = dist + 1
                        elif chain1[i] in ["A", "G"]:
                            dist = dist + 2
                    if chain2[i] == "S":
                        if chain1[i] in ["G", "C"]:
                            dist = dist + 1
                        elif chain1[i] in ["T", "A"]:
                            dist = dist + 2
                    if chain2[i] == "W":
                        if chain1[i] in ["A", "T"]:
                            dist = dist + 1
                        elif chain1[i] in ["C", "G"]:
                            dist = dist + 2
                    if chain2[i] == "K":
                        if chain1[i] in ["G", "T"]:
                            dist = dist + 1
                        elif chain1[i] in ["C", "A"]:
                            dist = dist + 2
                    if chain2[i] == "M":
                        if chain1[i] in ["A", "C"]:
                            dist = dist + 1
                        elif chain1[i] in ["T", "G"]:
                            dist = dist + 2
            elif chain1[i] in ["R", "Y", "S", "W", "K", "M"]:
                if chain2[i] in ["R", "Y", "S", "W", "K", "M"]:
                    if chain1[i] == chain2[i]:
                        dist = dist + 0
                    else:
                        dist = dist + 2
                elif chain2[i] in ["A", "C", "T", "G"]:
                    if chain1[i] == "R":
                        if chain2[i] in ["A", "G"]:
                            dist = dist + 1
                        elif chain2[i] in ["T", "C"]:
                            dist = dist + 2
                    if chain1[i] == "Y":
                        if chain2[i] in ["C", "T"]:
                            dist = dist + 1
                        elif chain2[i] in ["A", "G"]:
                            dist = dist + 2
                    if chain1[i] == "S":
                        if chain2[i] in ["G", "C"]:
                            dist = dist + 1
                        elif chain2[i] in ["A", "T"]:
                            dist = dist + 2
                    if chain1[i] == "W":
                        if chain2[i] in ["A", "T"]:
                            dist = dist + 1
                        elif chain2[i] in ["C", "G"]:
                            dist = dist + 2
                    if chain1[i] == "K":
                        if chain2[i] in ["G", "T"]:
                            dist = dist + 1
                        elif chain2[i] in ["C", "A"]:
                            dist = dist + 2
                    if chain1[i] == "M":
                        if chain2[i] in ["A", "C"]:
                            dist = dist + 1
                        if chain2[i] in ["G", "T"]:
                            dist = dist + 2
    return dist




# these are diff individuals
list1 = "XXXXRGCXYRXXYCYGXYYXXXYCRACXXZACRXYXYXXRGRYRYRXRXRYXXRXXYYXXXYRYXCXGMXZGRRXXCRRYXRYCYRRXCTXYXCCGASGYYGRWG"
list2 = 'XXXXGGYXTRXCCTTAXTCXXXYCRRTXXCXYAXYXYRXRARCRYRXRXGTXXRXRXYXXXYGCXAYAXXGRGRTXTGGYXGCCYARXTTAYXTYAGSGCTARWR'
# dist = 40

compare(list1, list2)
#40 #now 50

# these are the same individual in replicates
list1_30 = 'AXXRRRYXYGXYYYTRYTTXGXYYRRYYZYGYGGTYCRCRRRCRYRRRYAYCYGTRYYRXXYRYGAYRCRRRARYYYRXTXRCCYGRRYYRXXXXXXXXXXXXXX'
list2_02 = 'XXXXRRYXYGXXYYTRXTTXXXYYRGYXXYXTGXTXCXXXRRCRYXXRXAYXXGXRYYXXXYZYXAXRXXRRAXXXXRXZXRCCYGRXYXXCXTYAGSGCCAGWR'
# dist = 4

compare(list1_30,list2_02)
4 # now 2
# FABULOUS

# now let's use this function to compare every pair of sequences.
distance = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/dist_ind_id_dual.tsv"

individuals = []
for i in ind_gen_joined:
    individuals.append(i)

# as a dictionary
ind_genotypes_dual_dict = {}
for i in ind_gen_joined:
    ind_gens = []
    seq1 = ind_gen_joined[i]
    for j in ind_gen_joined:
        seq2 = ind_gen_joined[j]
        d = compare(seq1, seq2)
        ind_gens.append(str(d))
    ind_genotypes_dual_dict[i] = ind_gens

geno_dual = pd.DataFrame.from_dict(ind_genotypes_dual_dict)

# attempt to make a sequence for SIMP and one for LWED, since some SNPs are informative in one species and not the other.
#____________________________________******__________________________
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
                if (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) >= 7: # changed to 20 from 10
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
                elif (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3])) < 7: # changed from 10 to 20
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
ind_genotypes_dual_dict_lwed = {}
for i in ind_gen_joined_lwed:
    ind_gens = []
    seq1 = ind_gen_joined_lwed[i]
    for j in ind_gen_joined_lwed:
        seq2 = ind_gen_joined_lwed[j]
        d = compare(seq1, seq2)
        ind_gens.append(str(d))
    ind_genotypes_dual_dict_lwed[i] = ind_gens

geno_dual_lwed = pd.DataFrame.from_dict(ind_genotypes_dual_dict_lwed)
t_lwed= geno_dual_lwed.transpose().astype(int)
t_lwed.columns = individuals_lwed
new_t_lwed = t_lwed.rename(columns = {'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"}, index={'tamOpt-008_S8':"8", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':"27", 'tamOpt-026_S18':"26", 'tamOpt-031_S28':"31", 'tamOpt-002_S2':"2", 'tamOpt-010_S9':"10", 'tamOpt-003_S3':"3", \
    'tamOpt-033_S31':"33", 'tamOpt-012_S10':"12", 'tamOpt-030_S26':"30", 'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14"})


lwed_mask = np.zeros_like(new_t_lwed, dtype=np.bool)
lwed_mask[np.triu_indices_from(lwed_mask)] = True
lwed_mask[np.diag_indices_from(lwed_mask)] = False
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_lwed, mask=lwed_mask, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for L. weddelli - Dual and lwed SNP sets")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_LWED_duallwedsets_20threshold.png", dpi=300)

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
ind_genotypes_dual_dict_simp = {}
for i in ind_gen_joined_simp:
    ind_gens = []
    seq1 = ind_gen_joined_simp[i]
    for j in ind_gen_joined_simp:
        seq2 = ind_gen_joined_simp[j]
        d = compare(seq1, seq2)
        ind_gens.append(str(d))
    ind_genotypes_dual_dict_simp[i] = ind_gens

geno_dual_simp = pd.DataFrame.from_dict(ind_genotypes_dual_dict_simp)
t_simp= geno_dual_simp.transpose().astype(int)
t_simp.columns = individuals_simp

new_t_simp = t_simp.rename(columns = {'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"}, index={'tamOpt-001_S1':"1", 'tamOpt-016_S14':"16", 'tamOpt-018_S17':"18", 'tamOpt-017_S15':"17", \
       'tamOpt-020_S21':"20", 'tamOpt-029_S24':"29", 'tamOpt-023_S27':"23", 'tamOpt-015_S13':"15", \
       'tamOpt-006_S6':"6", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-024_S29':"24", \
       'tamOpt-021_S23':"21", 'tamOpt-028_S22':"28"})


simp_mask = np.zeros_like(new_t_simp, dtype=np.bool)
simp_mask[np.triu_indices_from(simp_mask)] = True
simp_mask[np.diag_indices_from(simp_mask)] = False
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(new_t_simp, mask=simp_mask, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for S. imperator - Dual and lwed SNP sets")
ax.tick_params(axis="y", labelrotation = 0)
plt.ylabel("Samples")
plt.xlabel("Samples")
ax.figure.tight_layout()






# okay so we have two new matrices. I should be able to get replicates, relatives and unique samples... right?
# let's start with SIMP
# the matrix is t_simp, and now we need the 
# Replicates in SIMP:
replicates_simp = ['tamOpt-001_S1','tamOpt-029_S24', 'tamOpt-006_S6', 'tamOpt-028_S22', 'tamOpt-034_S32']
rep_set_1_simp = ['tamOpt-001_S1','tamOpt-029_S24']
rep_set_2_simp = ['tamOpt-006_S6', 'tamOpt-028_S22', 'tamOpt-034_S32']
list_of_dist_replicates_simp = []
list_of_dist_replicates_simp.append(t_simp['tamOpt-001_S1']['tamOpt-029_S24'])
list_of_dist_replicates_simp.append(t_simp['tamOpt-006_S6']['tamOpt-028_S22'])
list_of_dist_replicates_simp.append(t_simp['tamOpt-006_S6']['tamOpt-034_S32'])
list_of_dist_replicates_simp.append(t_simp['tamOpt-028_S22']['tamOpt-034_S32'])

# Relatives in SIMP
unique_relatives_simp = ['tamOpt-018_S17', 'tamOpt-023_S27', 'tamOpt-020_S21', 'tamOpt-019_S19', 'tamOpt-017_S15', 'tamOpt-021_S23', 'tamOpt-024_S29']
all_relatives_simp = ['tamOpt-001_S1','tamOpt-029_S24', 'tamOpt-018_S17', 'tamOpt-023_S27', 'tamOpt-020_S21','tamOpt-019_S19', 'tamOpt-017_S15', 'tamOpt-006_S6', 'tamOpt-028_S22', 'tamOpt-034_S32','tamOpt-021_S23', 'tamOpt-024_S29' ]
rel_set_1_simp = ['tamOpt-018_S17', 'tamOpt-023_S27', 'tamOpt-020_S21']
rel_set_2_simp = ['tamOpt-019_S19', 'tamOpt-017_S15']
rel_set_3_simp = ['tamOpt-021_S23', 'tamOpt-024_S29']
list_of_dist_relatives_simp =[]
list_of_dist_relatives_simp.append(t_simp['tamOpt-023_S27']['tamOpt-018_S17'])
list_of_dist_relatives_simp.append(t_simp['tamOpt-020_S21']['tamOpt-018_S17'])
list_of_dist_relatives_simp.append(t_simp['tamOpt-020_S21']['tamOpt-023_S27'])
for i in rel_set_1_simp:
    for j in rep_set_1_simp:
        list_of_dist_relatives_simp.append(t_simp[i][j]) # end rel set 1

list_of_dist_relatives_simp.append(t_simp['tamOpt-019_S19']['tamOpt-017_S15'])
for i in rep_set_2_simp:
    for j in rel_set_2_simp:
        list_of_dist_relatives_simp.append(t_simp[i][j]) # end rel set 2

list_of_dist_relatives_simp.append(t_simp['tamOpt-021_S23']['tamOpt-024_S29']) # end rel set 3
# we expect 17 of these distances

# Unique samples in SIMP
list_of_dist_unique_simp = []
for i in rep_set_1_simp:
    for j in rep_set_2_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # replicate sets against eachother -IF NOT FROM SAME INDIVIDUAL

for i in rep_set_1_simp:
    for j in rel_set_2_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # replicate set 1 vs rel set 2

for i in rep_set_1_simp:
    for j in rel_set_3_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # replicate set 1 vs rel set 3

for i in rep_set_2_simp:
    for j in rel_set_1_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # replicate set 2 vs rel set 1

for i in rep_set_2_simp:
    for j in rel_set_3_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # replicate set 2 vs rel set 3

for i in rel_set_1_simp:
    for j in rel_set_2_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # rel set 1 vs rel set 2

for i in rel_set_1_simp:
    for j in rel_set_3_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # rel set 1 vs rel set 3

for i in rel_set_2_simp:
    for j in rel_set_3_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # rel set 2 vs rel set 3

unique_inds_simp = []
for ind in individuals_simp:
    if ind not in replicates_simp and ind not in all_relatives_simp:
        unique_inds_simp.append(ind)

list_of_dist_unique_simp.append(t_simp['tamOpt-016_S14']['tamOpt-015_S13']) # uniques vs uniques

for i in unique_inds_simp:
    for j in replicates_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # uniques vs all reps

for i in unique_inds_simp:
    for j in unique_relatives_simp:
        list_of_dist_unique_simp.append(t_simp[i][j]) # uniques vs all unique rels


# all values in the matrix are accounted for!!!!
# let's plot and test
rep_list = ["rep"]*4
unique_simp_list = ["unique"]*70
relative_list = ["relatives"]*17
status = rep_list + relative_list + unique_simp_list
distances = list_of_dist_replicates_simp + list_of_dist_relatives_simp + list_of_dist_unique_simp
distances_dict = {"status": status, "distance": distances}
df_simp = pd.DataFrame(distances_dict)
uniques_simp = df_simp.query("status == 'unique'")["distance"]
reps_simp = df_simp.query("status == 'rep'")["distance"]
rels_simp = df_simp.query("status == 'relatives'")["distance"]
df_simp.groupby("distance").describe()

# testing normality:
stats.shapiro(uniques_simp)
# ShapiroResult(statistic=0.9645323753356934, pvalue=0.04477005451917648) # not norm but close!
sns.histplot(uniques_simp)

stats.shapiro(rels)
# ShapiroResult(statistic=0.9388715028762817, pvalue=0.30504727363586426) # norm
stats.shapiro(reps)
# ShapiroResult(statistic=0.8949451446533203, pvalue=0.4063870310783386)

# equality of variances:
stats.levene(uniques_simp, reps, rels)
# LeveneResult(statistic=5.776584547835112, pvalue=0.004393530379604505) # not equal I think?

# welchs t test
def welch_ttest(x, y): 
    ## Welch-Satterthwaite Degrees of Freedom ##
    dof = (x.var()/x.size + y.var()/y.size)**2 / ((x.var()/x.size)**2 / (x.size-1) + (y.var()/y.size)**2 / (y.size-1))
    t, p = stats.ttest_ind(x, y, equal_var = False)
    print("\n",
          f"Welch's t-test= {t:.4f}", "\n",
          f"p-value = {p:.4f}", "\n",
          f"Welch-Satterthwaite Degrees of Freedom= {dof:.4f}")

welch_ttest(uniques_simp,reps)
#  Welch's t-test= 30.7349 
# p-value = 0.0000 
# Welch-Satterthwaite Degrees of Freedom= 38.6896
welch_ttest(uniques_simp,rels)
# Welch's t-test= 4.2771 
# p-value = 0.0004 
# Welch-Satterthwaite Degrees of Freedom= 19.9767

from statannot import add_stat_annotation
ax = sns.boxplot(x="status", y = "distance", data = df_simp)
add_stat_annotation(ax, data = df_simp, x = "status", y = "distance",  box_pairs = [("unique", "rep"), ("unique", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/rep_vs_rel_vs_unique_SIMP.png", dpi=300)

# LWED
# Replicates in LWED
under_performing_inds = ["tamOpt-035_S33", "tamOpt-007_S7", 'tamOpt-022_S25', 'tamOpt-025_S16', "tamOpt-004_S4"]
replicates_lwed = ['tamOpt-002_S2', 'tamOpt-030_S26', 'tamOpt-003_S3', 'tamOpt-031_S28','tamOpt-026_S18', 'tamOpt-032_S30', 'tamOpt-005_S5', 'tamOpt-027_S20','tamOpt-033_S31']
rep_set_1_lwed = ['tamOpt-002_S2', 'tamOpt-030_S26']
rep_set_2_lwed = ['tamOpt-003_S3', 'tamOpt-031_S28']
rep_set_3_lwed = ['tamOpt-026_S18', 'tamOpt-032_S30']
rep_set_4_lwed = ['tamOpt-005_S5', 'tamOpt-027_S20', 'tamOpt-033_S31']
list_of_dist_replicates_lwed = []
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-002_S2']['tamOpt-030_S26'])
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-003_S3']['tamOpt-031_S28'])
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-026_S18']['tamOpt-032_S30'])
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-005_S5']['tamOpt-027_S20'])
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-005_S5']['tamOpt-033_S31'])
list_of_dist_replicates_lwed.append(t_lwed['tamOpt-027_S20']['tamOpt-033_S31'])
# we got 6 of these

# Relatives in LWED
unique_relatives_lwed = []
all_relatives_lwed = ['tamOpt-003_S3','tamOpt-031_S28','tamOpt-005_S5', 'tamOpt-027_S20', 'tamOpt-033_S31']
rel_set_1_lwed = [] # look above, they're all related 
list_of_dist_relatives_lwed = []
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-003_S3']['tamOpt-005_S5'])
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-003_S3']['tamOpt-027_S20'])
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-003_S3']['tamOpt-033_S31'])
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-031_S28']['tamOpt-005_S5'])
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-031_S28']['tamOpt-027_S20'])
list_of_dist_relatives_lwed.append(t_lwed['tamOpt-031_S28']['tamOpt-033_S31'])

# Uniques in LWED
unique_inds_lwed = []
for ind in individuals_lwed:
    if ind not in replicates_lwed and ind not in all_relatives_lwed and ind not in under_performing_inds:
        unique_inds_lwed.append(ind)

list_of_dist_unique_lwed = []
for i in unique_inds_lwed:
    for j in replicates_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # uniques vs replicates

for i in rep_set_1_lwed:
    for j in rep_set_2_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # rep set 1 vs rep set 2

for i in rep_set_1_lwed:
    for j in rep_set_3_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # rep set 1 vs rep set 3

for i in rep_set_1_lwed:
    for j in rep_set_4_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # rep set 1 vs rep set 4

for i in rep_set_2_lwed:
    for j in rep_set_3_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # rep set 2 vs rep set 3

for i in rep_set_3_lwed:
    for j in rep_set_4_lwed:
        list_of_dist_unique_lwed.append(t_lwed[i][j]) # rep set 3 vs rep set 4

from itertools import combinations
for comb in combinations(unique_inds_lwed, 2):
    ind1 = comb[0]
    ind2 = comb[1]
    list_of_dist_unique_lwed.append(t_lwed[ind1][ind2]) # uniques vs uniques

# exactly the number of distance values we expected! 105
# let's plot and test
rep_list = ["rep"]*6
unique_lwed_list = ["unique"]*79
relative_list = ["relatives"]*6
status = rep_list + relative_list + unique_lwed_list
distances = list_of_dist_replicates_lwed + list_of_dist_relatives_lwed + list_of_dist_unique_lwed
distances_dict = {"status": status, "distance": distances}
df_lwed = pd.DataFrame(distances_dict)
uniques_lwed = df_lwed.query("status == 'unique'")["distance"]
reps_lwed = df_lwed.query("status == 'rep'")["distance"]
rels_lwed = df_lwed.query("status == 'relatives'")["distance"]
df_lwed.groupby("distance").describe()

# testing normality:
stats.shapiro(uniques_lwed)
# ShapiroResult(statistic=0.9886779189109802, pvalue=0.7172327041625977) # norm 
sns.histplot(uniques_lwed)

stats.shapiro(rels)
# ShapiroResult(statistic=0.859739363193512, pvalue=0.1882738471031189) # norm
stats.shapiro(reps)
#ShapiroResult(statistic=0.9024807214736938, pvalue=0.3888012170791626) # norm

# equality of variances:
stats.levene(uniques_lwed, reps, rels)
# LeveneResult(statistic=4.569506499072931, pvalue=0.012939436757058033) # not equal I think?

welch_ttest(uniques_lwed,reps)
# Welch's t-test= 25.6589 
# p-value = 0.0000 
# Welch-Satterthwaite Degrees of Freedom= 36.8434

welch_ttest(uniques_lwed,reps)
# Welch's t-test= 25.6589 
# p-value = 0.0000 
# Welch-Satterthwaite Degrees of Freedom= 36.8434

from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


ax = sns.boxplot(x="status", y = "distance", data = df_lwed)
add_stat_annotation(ax, data = df_lwed, x = "status", y = "distance",  box_pairs = [("unique", "rep"), ("unique", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/rep_vs_rel_vs_unique_LWED.png", dpi=300)

# plot both df_simp and df_lwed sidebyside:
fig, (ax1, ax2) = plt.subplots(ncols=2)
sns.boxplot(x="status", y = "distance", data = df_lwed, ax=ax1)
ax1.set_ylabel("Distance")
ax1.set_xlabel("L. weddelli samples")
add_stat_annotation(ax1, data = df_lwed, x = "status", y = "distance",  box_pairs = [("unique", "rep"), ("unique", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)

sns.boxplot(x="status", y = "distance", data = df_simp, ax = ax2)
ax2.set_xlabel("S. imperator samples")
add_stat_annotation(ax2, data = df_simp, x = "status", y = "distance",  box_pairs = [("unique", "rep"), ("unique", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)
ax2.set(ylim=(-4,119))
ax2.set(yticklabels=[])
ax2.set_ylabel("")
ax1.figure.tight_layout()
ax2.figure.tight_layout()




plt.subplots_adjust(wspace=0.1)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/rep_vs_rel_vs_unique_both_separate_report.png", dpi=300)








# okay so let's see if we can add the reps and rels and just plot it all together?

rep_list = ["rep"]*10
unique_lwed_list = ["unique"]*149
relative_list = ["relatives"]*23
status = rep_list + relative_list + unique_lwed_list
distances = list_of_dist_replicates_lwed + list_of_dist_replicates_simp + list_of_dist_relatives_lwed + list_of_dist_relatives_simp + list_of_dist_unique_lwed + list_of_dist_unique_simp
distances_dict = {"status": status, "distance": distances}
df = pd.DataFrame(distances_dict)
uniques = df.query("status == 'unique'")["distance"]
reps = df.query("status == 'rep'")["distance"]
rels = df.query("status == 'relatives'")["distance"]
df.groupby("distance").describe()
import scipy.stats as stats

stats.shapiro(uniques)
# ShapiroResult(statistic=0.9929917454719543, pvalue=0.6827273368835449) norm
stats.shapiro(reps)
# ShapiroResult(statistic=0.8902443051338196, pvalue=0.17065417766571045) norm
stats.shapiro(rels)
# ShapiroResult(statistic=0.9167919158935547, pvalue=0.056893669068813324) just at the edge! but norm
stats.levene(uniques, reps, rels)
# LeveneResult(statistic=9.00792662588885, pvalue=0.00018727991752333112) # not equal so we have to do welchs
from statannot import add_stat_annotation
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
sns.set(font="serif")
ax = sns.boxplot(x="status", y = "distance", data = df)
add_stat_annotation(ax, data = df, x = "status", y = "distance",  box_pairs = [("unique", "rep"), ("unique", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/rep_vs_rel_vs_unique_BOTH.png", dpi=300, transparent= True)








#____________________________________******__________________________


# let's make one per species _________________________
under_performing_inds = ["tamOpt-035_S33", "tamOpt-007_S7", 'tamOpt-022_S25', 'tamOpt-025_S16', "tamOpt-004_S4"]

# SIMP
ind_genotypes_dual_dict_SIMP = {}
for i in ind_gen_joined:
    ind_gens = []
    seq1 = ind_gen_joined[i]
    if ind_spp[i] == "SIMP" and i not in under_performing_inds:
        for j in ind_gen_joined:
            if ind_spp[j] == "SIMP" and j not in under_performing_inds:
                seq2 = ind_gen_joined[j]
                d = compare(seq1, seq2)
                ind_gens.append(str(d))
        ind_genotypes_dual_dict_SIMP[i] = ind_gens

geno_dual_SIMP = pd.DataFrame.from_dict(ind_genotypes_dual_dict_SIMP)
inds = []
for col in geno_dual_SIMP.columns:
    inds.append(col)

t_simp = geno_dual_SIMP.transpose().astype(int)
t_simp.columns = inds

# LWED
ind_genotypes_dual_dict_LWED = {}
for i in ind_gen_joined:
    ind_gens = []
    seq1 = ind_gen_joined[i]
    if ind_spp[i] == "LWED" and i not in under_performing_inds:
        for j in ind_gen_joined:
            if ind_spp[j] == "LWED" and j not in under_performing_inds:
                seq2 = ind_gen_joined[j]
                d = compare(seq1, seq2)
                ind_gens.append(str(d))
        ind_genotypes_dual_dict_LWED[i] = ind_gens

geno_dual_LWED = pd.DataFrame.from_dict(ind_genotypes_dual_dict_LWED)
inds = []
for col in geno_dual_LWED.columns:
    inds.append(col)

t_lwed = geno_dual_LWED.transpose().astype(int)
t_lwed.columns = inds

sns.clustermap(t_lwed, annot= True)

# now let's try and plot these nicely
simp_mask = np.zeros_like(t_simp, dtype=np.bool)
simp_mask[np.triu_indices_from(simp_mask)] = True
simp_mask[np.diag_indices_from(simp_mask)] = False
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(t_simp, mask=simp_mask, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for S. imperator")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_SIMP.png", dpi=300)

lwed_mask = np.zeros_like(t_lwed, dtype=np.bool)
lwed_mask[np.triu_indices_from(lwed_mask)] = True
lwed_mask[np.diag_indices_from(lwed_mask)] = False
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(t_lwed, mask=lwed_mask, annot= True, annot_kws={"size":10}, cmap="Blues")
ax.set_title("Distance matrix of custom SNP-only sequences for L. weddelli")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/distance_heatmap_LWED.png", dpi=300)



#___________________________________________
t= geno_dual.transpose().astype(int)
t.columns = individuals
p = sns.clustermap(t, annot= T)

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
plt.show()

import numpy as np

# let's try and annotate the spp and ind id of these.
spp_id = ['SIMP', 'SIMP', 'LWED', 'SIMP', 'NA', 'LWED', 'LWED', 'LWED', 'SIMP', 'LWED', 'LWED', 'SIMP', 'NA', 'SIMP', 'LWED', 'SIMP', 'SIMP', 'LWED', 'LWED', 'SIMP', 'LWED', 'SIMP', 'SIMP', 'LWED', 'SIMP', 'NA', 'NA', 'SIMP', 'LWED', 'LWED', 'LWED', 'SIMP', 'LWED']
# this didn't work: trying to do the same as naming clusters but for annotating samples with the individual so we could tell repeat samples spp_id = ["73", "644", "326", "466", "BL", "10", "27", "10", "27", "10", "516", "216", "203", "1", "10", "73", "383", "1223", "610", "216", "27", "182", "388", "429", "182", "203", "159", "540", "223", "524", "75", "27", "148", "182", "223"]
net_spp_id = sorted(np.unique(spp_id, return_index=True)[1])
spp_names =[spp_id[index] for index in sorted(net_spp_id)]
spp_pal = sns.husl_palette(len(spp_names))
spp_lut = dict(zip(map(str, np.unique(spp_id)), spp_pal))
spp_colors = pd.Series(spp_id).map(spp_lut)
spp_colors = np.asarray(spp_colors)

ind_names = individuals
n_ind = len(ind_names)
M = np.random.rand(n_ind, n_ind)

p = sns.clustermap(t, row_colors = spp_colors, col_colors=spp_colors)
#p.ax_heatmap.set_xticklabels(ind_names, rotation = 90)
#p.ax_heatmap.set_yticklabels(ind_names, rotation = 0)
#p.ax_row_colors.set_yticks(0.5 * (np.array(net_spp_id) + np.array(net_spp_id[1:] + [len(spp_id)])))

# let's see if we can notice any differences in distance values for replicates vs unique samples.
replicates = ["tamOpt-001_S1", "tamOpt-029_S24", "tamOpt-002_S2", "tamOpt-030_S26", "tamOpt-003_S3", "tamOpt-031_S28", "tamOpt-026_S18", "tamOpt-032_S30", "tamOpt-005_S5", "tamOpt-027_S20", "tamOpt-033_S31", "tamOpt-006_S6", "tamOpt-028_S22", "tamOpt-034_S32"]
replicates_simp = ["tamOpt-001_S1", "tamOpt-029_S24","tamOpt-006_S6", "tamOpt-028_S22", "tamOpt-034_S32" ]
replicates_lwed = ["tamOpt-002_S2", "tamOpt-030_S26", "tamOpt-026_S18", "tamOpt-032_S30", "tamOpt-005_S5", "tamOpt-027_S20", "tamOpt-033_S31", "tamOpt-003_S3", "tamOpt-031_S28"]
rep_set_1 = ["tamOpt-001_S1", "tamOpt-029_S24"]
rep_set_2 = ["tamOpt-002_S2", "tamOpt-030_S26"]
rep_set_3 = ["tamOpt-003_S3", "tamOpt-031_S28"]
rep_set_4 = ["tamOpt-026_S18", "tamOpt-032_S30"]
rep_set_5 = ["tamOpt-005_S5", "tamOpt-027_S20", "tamOpt-033_S31"]
rep_set_6 = ["tamOpt-006_S6", "tamOpt-028_S22", "tamOpt-034_S32"]

# remove failed inds.
under_performing_inds = ["tamOpt-035_S33", "tamOpt-007_S7", 'tamOpt-022_S25', 'tamOpt-025_S16', "tamOpt-004_S4"]

# REPLICATES
list_of_dist_replicates = []

list_of_dist_replicates.append(t['tamOpt-001_S1']['tamOpt-029_S24'])
list_of_dist_replicates.append(t['tamOpt-030_S26']['tamOpt-002_S2'])
list_of_dist_replicates.append(t['tamOpt-003_S3']['tamOpt-031_S28'])
list_of_dist_replicates.append(t['tamOpt-032_S30']['tamOpt-026_S18'])
list_of_dist_replicates.append(t['tamOpt-027_S20']['tamOpt-033_S31'])
list_of_dist_replicates.append(t['tamOpt-027_S20']['tamOpt-005_S5'])
list_of_dist_replicates.append(t['tamOpt-033_S31']['tamOpt-005_S5'])
list_of_dist_replicates.append(t['tamOpt-034_S32']['tamOpt-028_S22'])
list_of_dist_replicates.append(t['tamOpt-034_S32']['tamOpt-006_S6'])
list_of_dist_replicates.append(t['tamOpt-006_S6']['tamOpt-028_S22'])


# RELATIVES (full siblings or parent-offspring pairs)
unique_relatives = ['tamOpt-023_S27', 'tamOpt-020_S21','tamOpt-019_S19', 'tamOpt-017_S15','tamOpt-021_S23','tamOpt-024_S29', 'tamOpt-018_S17']
rel_set_1 = ['tamOpt-018_S17','tamOpt-023_S27', 'tamOpt-020_S21'] #rep_set_1[0], rep_set_1[1]] 
rel_set_2 = rep_set_3 + rep_set_5 
rel_set_3 = ['tamOpt-019_S19', 'tamOpt-017_S15'] #rep_set_6[0], rep_set_6[1], rep_set_6[2]]
rel_set_4 = ['tamOpt-021_S23','tamOpt-024_S29']
relatives_simp = ['tamOpt-023_S27', 'tamOpt-020_S21', 'tamOpt-019_S19', 'tamOpt-017_S15','tamOpt-021_S23','tamOpt-024_S29']
relatives_lwed = []


list_of_dist_relatives =[]
list_of_dist_relatives.append(t['tamOpt-023_S27']['tamOpt-018_S17'])
list_of_dist_relatives.append(t['tamOpt-001_S1']['tamOpt-018_S17'])
list_of_dist_relatives.append(t['tamOpt-029_S24']['tamOpt-018_S17'])
list_of_dist_relatives.append(t['tamOpt-020_S21']['tamOpt-018_S17'])
list_of_dist_relatives.append(t['tamOpt-023_S27']['tamOpt-020_S21'])
list_of_dist_relatives.append(t['tamOpt-023_S27']['tamOpt-001_S1'])
list_of_dist_relatives.append(t['tamOpt-023_S27']['tamOpt-029_S24'])
list_of_dist_relatives.append(t['tamOpt-020_S21']['tamOpt-001_S1'])
list_of_dist_relatives.append(t['tamOpt-020_S21']['tamOpt-029_S24']) # end rel_set1_1

for i in rep_set_3:
    for j in rep_set_5:
        list_of_dist_relatives.append(t[i][j])  # end rel_set_2

for i in rel_set_3:
    for j in rep_set_6:
        list_of_dist_relatives.append(t[i][j]) # end rel_set_3

list_of_dist_relatives.append(t['tamOpt-021_S23']['tamOpt-024_S29']) # end_rel_set_4

# list is 22 values long, EXACTLY what we expected.


## UNIQUE
replicates = ["tamOpt-001_S1", "tamOpt-029_S24", "tamOpt-002_S2", "tamOpt-030_S26", "tamOpt-003_S3", "tamOpt-031_S28", "tamOpt-026_S18", "tamOpt-032_S30", "tamOpt-005_S5", "tamOpt-027_S20", "tamOpt-033_S31", "tamOpt-006_S6", "tamOpt-028_S22", "tamOpt-034_S32"]
unique_relatives = ['tamOpt-023_S27', 'tamOpt-020_S21','tamOpt-019_S19', 'tamOpt-017_S15','tamOpt-021_S23','tamOpt-024_S29', 'tamOpt-018_S17']
under_performing_inds = ["tamOpt-035_S33", "tamOpt-007_S7", 'tamOpt-022_S25', 'tamOpt-025_S16', "tamOpt-004_S4"]
unique = []
unique_simp = []
unique_lwed = []
for ind in individuals:
    if ind_spp[ind] == "SIMP":
        if ind not in replicates and ind not in unique_relatives and ind not in under_performing_inds:
            unique_simp.append(ind)
    elif ind_spp[ind] == "LWED":
        if ind not in replicates and ind not in unique_relatives and ind not in under_performing_inds:
            unique_lwed.append(ind)

from itertools import combinations
list_of_dist_unique_simp =[]
for comb in combinations(unique_simp, 2):
    ind1 = comb[0]
    ind2 = comb[1]
    list_of_dist_unique_simp.append(t[ind1][ind2])

list_of_dist_unique_lwed =[]
for comb in combinations(unique_lwed, 2):
    ind1 = comb[0]
    ind2 = comb[1]
    list_of_dist_unique_lwed.append(t[ind1][ind2])

# plus the unique samples with the replicates:

for i in unique_simp:
    for j in replicates_simp:
        list_of_dist_unique_simp.append(t[i][j])

for i in unique_lwed:
    for j in replicates_lwed:
        list_of_dist_unique_lwed.append(t[i][j])


for i in unique_simp:
    for j in relatives_simp:
        list_of_dist_unique_simp.append(t[i][j])

for i in unique_lwed:
    for j in relatives_lwed:
        list_of_dist_unique_lwed.append(t[i][j])


for i in replicates_simp:
    for j in relatives_simp:
        list_of_dist_unique_simp.append(t[i][j])

for i in replicates_lwed:
    for j in relatives_lwed:
        list_of_dist_unique_lwed.append(t[i][j])


# ready!



rep_list = ["rep"]*10

#unique_list = ["unique"]*272
unique_simp_list = ["unique_simp"]*53
unique_lwed_list = ["unique_lwed"]*55

relative_list = ["relatives"]*22

status = rep_list + relative_list + unique_simp_list + unique_lwed_list

distances = list_of_dist_replicates + list_of_dist_relatives + list_of_dist_unique_simp + list_of_dist_unique_lwed

distances_dict = {"status": status, "distance": distances}

df = pd.DataFrame(distances_dict)

df.boxplot(by = "status", column = ["distance"])



import scipy.stats as stats

stats.ttest_ind(df["distance"][df["status"] == "rep"], df["distance"][df["status"] == "unique_simp"])
#Ttest_indResult(statistic=-8.155314779811611, pvalue=2.378255833217758e-11)
stats.ttest_ind(df["distance"][df["status"] == "rep"], df["distance"][df["status"] == "unique_lwed"])
#Ttest_indResult(statistic=-10.764993827217346, pvalue=6.522282433274143e-16)
stats.ttest_ind(df["distance"][df["status"] == "rep"], df["distance"][df["status"] == "relatives"])
#Ttest_indResult(statistic=-6.406030348091335, pvalue=4.508199830083702e-07)
stats.ttest_ind(df["distance"][df["status"] == "relatives"], df["distance"][df["status"] == "unique_simp"])
#Ttest_indResult(statistic=-1.9029126313857498, pvalue=0.06099557830649008)
stats.ttest_ind(df["distance"][df["status"] == "relatives"], df["distance"][df["status"] == "unique_lwed"])
#Ttest_indResult(statistic=-4.326285060359399, pvalue=4.614611338142373e-05)


ax = sns.boxplot(data=df, y = "distance", x = "status")

#uniques = df.query("status == 'unique'")["distance"]
uniques_simp = df.query("status == 'unique_simp'")["distance"]
uniques_lwed = df.query("status == 'unique_lwed'")["distance"]
reps = df.query("status == 'rep'")["distance"]
rels = df.query("status == 'relatives'")["distance"]
df.groupby("distance").describe()

# check normality of the three sets
stats.shapiro(uniques_simp)
#ShapiroResult(statistic=0.9636923670768738, pvalue=0.10685442388057709) # normally distributed
sns.histplot(uniques_simp)

stats.shapiro(uniques_lwed)
#ShapiroResult(statistic=0.9865913987159729, pvalue=0.7959246039390564) # normally distributed
sns.histplot(uniques_lwed)

stats.shapiro(reps)
#ShapiroResult(statistic=0.8902443051338196, pvalue=0.17065417766571045) # normally distributed
sns.histplot(reps)

stats.shapiro(rels)
#ShapiroResult(statistic=0.9185061454772949, pvalue=0.07079616189002991) # normally distributed
sns.histplot(rels)


# make sure variances are equal
stats.levene(uniques_simp, reps, rels)
#LeveneResult(statistic=7.7349544685974125, pvalue=0.0008369027109971675) # no
stats.levene(uniques_lwed, reps, rels)
#LeveneResult(statistic=7.471458743919302, pvalue=0.0010317443694725058) # no
stats.levene(uniques_lwed, uniques_simp)
#LeveneResult(statistic=0.813395686803928, pvalue=0.36916295720475234) # equal


import statistics
statistics.variance(uniques_simp)
#184.99419448476053
statistics.variance(uniques_lwed)
#155.84309764309765
statistics.variance(reps)
#4.011111111111111
statistics.variance(rels)
#196.16450216450215





# t- test assuming we meet the requirements
res = stats.ttest_ind(uniques,reps, equal_var = True)


# welchs t test if variances are not equal:
res = stats.ttest_ind(uniques_simp,reps, equal_var = False)
#Ttest_indResult(statistic=17.932463246619626, pvalue=8.545716454702923e-26)
# let's try and get the degrees of freedom:
def welch_dof(x,y):
    dof = (x.var()/x.size + y.var()/y.size)**2 / ((x.var()/x.size)**2 / (x.size-1) + (y.var()/y.size)**2 / (y.size-1))
    print(f"Welch-Satterthwaite Degrees of Freedom= {dof:.4f}")

welch_dof(uniques_simp,reps)
# Welch-Satterthwaite Degrees of Freedom= 73.8665
def welch_ttest(x, y): 
    ## Welch-Satterthwaite Degrees of Freedom ##
    dof = (x.var()/x.size + y.var()/y.size)**2 / ((x.var()/x.size)**2 / (x.size-1) + (y.var()/y.size)**2 / (y.size-1))
    t, p = stats.ttest_ind(x, y, equal_var = False)
    print("\n",
          f"Welch's t-test= {t:.4f}", "\n",
          f"p-value = {p:.4f}", "\n",
          f"Welch-Satterthwaite Degrees of Freedom= {dof:.4f}")

welch_ttest(uniques_simp,reps)
# Welch's t-test= 17.9325 
# p-value = 0.0000 
# Welch-Satterthwaite Degrees of Freedom= 60.0558

welch_ttest(uniques_lwed,rels)
# Welch's t-test= 4.1160 
# p-value = 0.0002 
# Welch-Satterthwaite Degrees of Freedom= 35.0894

welch_ttest(uniques_simp,rels)
# Welch's t-test= 1.8797 
# p-value = 0.0678 
# Welch-Satterthwaite Degrees of Freedom= 38.2897

welch_ttest(uniques_lwed,rels)
# Welch's t-test= 4.1160 
# p-value = 0.0002 
# Welch-Satterthwaite Degrees of Freedom= 35.0894

from scipy.stats import f_oneway
f_oneway(uniques_lwed,uniques_simp, rels, reps)
# F_onewayResult(statistic=33.72377080740251, pvalue=2.326378671037741e-16)

import matplotlib.pyplot as plt
import pandas as pd
from utils import *
import numpy as np
from scipy.stats import mannwhitneyu, normaltest


from statsmodels.stats.weightstats import ttest_ind
ttest_ind(uniques_simp, reps)
# (8.15531477981161, 2.3782558332177747e-11, 61.0)
ttest_ind(uniques_simp, rels)


statistics.stdev(uniques) #  SD = 24.554452209754615  # mean = 69.96703296703296
statistics.stdev(reps) # SD = 1.9021518914592017  #  mean =  2.727272727272727
statistics.mean(uniques)


from statannot import add_stat_annotation

ax = sns.boxplot(x="status", y = "distance", data = df)
add_stat_annotation(ax, data = df, x = "status", y = "distance",  box_pairs = [("unique_simp", "rep"), ("unique_lwed", "rep"), ("unique_simp", "relatives"), ("unique_lwed", "relatives"), ("relatives", "rep")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plt_unique_rep.png", dpi=300)

plt.show()

