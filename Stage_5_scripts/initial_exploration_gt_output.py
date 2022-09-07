import csv
import glob
from os import spawnlpe
from traceback import format_exception_only
import pandas as pd
from pandas import DataFrame
from collections import defaultdict

print(glob.glob("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/updated_geno/*.genos"))

dict_of_inds = {}

# for all output files I forgot to update the first_geno dir to second_geno run, will have to rerun everything later

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

# okay, we now have a dictionary containing all individuals, in each individual is another dictionary for each SNP (using the bp3 index), the components of that
# dictionary is a list of the a1, the a1counts, the a2, the a2counts, the genotype, and the genotype class.

# in case I need the snp dir being the main one
snp_dir = DataFrame(dict_of_inds).transpose().to_dict()


# okay, so firstly, I'd like to see how many times each SNP failed to get sequenced.

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

failed_ind_id = []
failed_lwed = []
failed_simp = []
failed_sex = []
failed_spp = []

for snp in no_reads_snps:
    if snp in ind_id_dual:
        failed_ind_id.append(snp)
    if snp in ind_id_lwed:
        failed_lwed.append(snp)
    if snp in ind_id_simp:
        failed_simp.append(snp)
    if snp in sex_id:
        failed_sex.append(snp)
    if snp in species_id:
        failed_spp.append(snp)

# produced reads (regardless of how many)

prod_reads= []

for snp in snp_dir:
    if snp not in no_reads_snps:
        prod_reads.append(snp)

prod_reads_dual = []
prod_reads_lwed = []
prod_reads_simp= []
prod_reads_spp = []
prod_reads_sex = []

for snp in prod_reads:
    if snp in ind_id_dual:
        prod_reads_dual.append(snp)
    if snp in ind_id_lwed:
        prod_reads_lwed.append(snp)
    if snp in ind_id_simp:
        prod_reads_simp.append(snp)
    if snp in sex_id:
        prod_reads_sex.append(snp)
    if snp in species_id:
        prod_reads_spp.append(snp)

# bp3 error

from Bio import SeqIO
import re
primers_file="/proj/proj_name/nobackup/SAM/flanking_regions/final_merged_fastas/19.02.2022_fastas_for_merging/simulating_seqs_to_troubleshoot_gt_seq/primers.tsv"
fasta_new_snps = "/proj/proj_name/nobackup/SAM/flanking_regions/final_merged_fastas/19.02.2022_fastas_for_merging/simulating_seqs_to_troubleshoot_gt_seq/concatenated_ALL_19.02.2022_adding_brakets_replace_snps.fasta"
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


# replacement SNPs with 33 genotyped inds, now these are 40 SNPs


# let's get stats for coverage for each SNP (avg cov, 1st quartile etc)
snp_dir_with_cov = {}
for snp in snp_dir:
    ind_dir = {}
    for ind in snp_dir[snp]:
        tuple_to_list = list(snp_dir[snp][ind])
        total_cov = int(tuple_to_list[1]) + int(tuple_to_list[3])
        if ind not in ind_dir:
            ind_dir[ind] = snp_dir[snp][ind][0], snp_dir[snp][ind][1], snp_dir[snp][ind][2], snp_dir[snp][ind][3], snp_dir[snp][ind][4], snp_dir[snp][ind][5], total_cov
    if snp not in snp_dir_with_cov:
        snp_dir_with_cov[snp] = ind_dir

#creating a list of the individuals
ind_list=[]
for snp in snp_dir:
    for ind in snp_dir[snp]:
        ind_list.append(ind)
    break

# we now have a dictionary that contains the total coverage of each snp for each individual. 
dict_cov_per_snp = {}
for snp in snp_dir_with_cov:
    cov_list_inds = []
    for ind in snp_dir_with_cov[snp]:
        cov_list_inds.append(snp_dir_with_cov[snp][ind][-1])
    #if str(snp) not in snps_that_arent_working and snp not in no_reads_snps:
        dict_cov_per_snp[snp] = cov_list_inds


# okay I'd like to visualize a boxplot of each snp. I will make a df of the coverage for each snp.
import matplotlib.pyplot as plt
import numpy as np
#snp_cov_df = pd.DataFrame(dict_cov_per_snp)
#plt.boxplot(snp_cov_df)
#plt.show()

# even after removing the 66 SNPs that didn't generate any reads, it's still a bit messy to visualize all together. What I will do is generate three dataframes,
# one per each snp category. For this, I need 3 lists that contain the BP3 indexes for the sets.
ind_id_dual = [21,22,23,24,25,26,29,30,31,32,34,37,38,39,42,43,44,46,47,49,51,53,54,55,58,60,61,62,64,65,67,69,70,74,75,76,78,80,81,82,83,86,88,89,90,91,92,94,96,98,100,101,103,104,105,107,108,109,110,111,112,113,114,117,118,119,120,123,124,128,129,130,131,132,133,134,135,137,138,139,140,141,142,144,145,146,147,148,149,150,152,153,154,156,157,158,159,160,161,162,164,165,166,167,168,169,171,173,174,176,177,178,179,181,183,184,187,188,189,190,192,193,358,360,361,363,364,367,369,372,373,376,378,379,380,381,382,383,384,385,388,389,390,392,393,394,395,396,399,400,401]
ind_id_lwed = [226,227,228,230,231,232,233,234,235,236,237,239,240,241,242,243,244,245,246,248,249,250,251,252,253,254,255,256,257,258,259,260,261,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291]
ind_id_simp = [292,293,295,296,297,299,300,302,303,304,306,307,308,310,311,312,313,315,316,317,318,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,341,343,344,345,346,347,348,349,352,353,355,356,357]
sex_id = [195,196,197,198,200,201,203,205,208,209,211,213,215,216,218,219,220,222]
species_id = [1,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

ind_id_dir_cov = {}
ind_id_lwed_cov = {}
ind_id_simp_cov = {}
sex_id_cov = {}
species_id_cov = {}

for snp in dict_cov_per_snp:
    if snp in ind_id_dual:
        ind_id_dir_cov[snp] = dict_cov_per_snp[snp]
    if snp in ind_id_lwed:
        ind_id_lwed_cov[snp] = dict_cov_per_snp[snp]
    if snp in ind_id_simp:
        ind_id_simp_cov[snp] = dict_cov_per_snp[snp]
    if snp in sex_id:
        sex_id_cov[snp] = dict_cov_per_snp[snp]
    if snp in species_id:
        species_id_cov[snp] = dict_cov_per_snp[snp]

# DUAL SPECIES INDIVIDUAL IDENTIFICATION SNPS

dual_ind = pd.DataFrame.from_dict(ind_id_dir_cov, orient='columns')
dual_ind_1 = dual_ind.reindex(dual_ind.mean().sort_values().index,axis = 1)
list_snps = []
for key in dual_ind_1:
    list_snps.append(str(key))

list_snps = []
for key in dual_ind_1.iloc[:, -22:]:
    list_snps.append(str(key))


plt.boxplot(dual_ind_1.iloc[:, -22:])
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=10)
plt.show()

# LWED INDIVIDUAL IDENTIFICATION SNPs
lwed_ind = pd.DataFrame.from_dict(ind_id_lwed_cov)
lwed_ind_1 = lwed_ind.reindex(lwed_ind.mean().sort_values().index,axis = 1)
list_snps = []
for key in lwed_ind_1:
    list_snps.append(str(key))


plt.boxplot(lwed_ind_1)
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=6)
plt.show()

# SIMP INDIVIDUAL IDENTIFICATION SNPs
simp_ind = pd.DataFrame.from_dict(ind_id_simp_cov)
simp_ind_1 = simp_ind.reindex(simp_ind.mean().sort_values().index,axis = 1)

list_snps = []
for key in simp_ind_1:
    list_snps.append(str(key))


plt.boxplot(simp_ind_1)
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=10)
plt.show()


# SEX IDENTIFICATION SNPS
sex = pd.DataFrame.from_dict(sex_id_cov)
sex_1 = sex.reindex(sex.mean().sort_values().index,axis = 1)
list_snps = []
for key in sex_1:
    list_snps.append(str(key))


plt.boxplot(sex_1)
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=10)
plt.show()

# SPECIES IDENTIFICATION SNPS
species = pd.DataFrame.from_dict(species_id_cov)
species_1 = species.reindex(species.mean().sort_values().index,axis = 1)
list_snps = []
for key in species_1:
    list_snps.append(str(key))


plt.boxplot(species_1)
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=10)
plt.show()

# OKAY so for the next step I want to perhaps highlight the ones that failed at the primer design step. Perhaps they amplified but they are not informative.
# we have those in a list: snps_that_arent_working

for snp in dual_ind_1:
    if str(snp) in snps_that_arent_working:
        print(snp)

# dual ind = 390, 393 
# NOW only 390

lwed_ind_list = []
for snp in lwed_ind_1:
    if str(snp) in snps_that_arent_working:
        lwed_ind_list.append(snp)

#>>> lwed_ind_list
#[280, 269, 242, 236, 251, 265, 248, 261, 250, 244, 267, 272, 243, 230, 253, 241, 252, 234, 260, 291, 266, \
#    275, 256, 258, 283, 279, 245, 273, 228, 240, 263, 246, 255, 227, 231, 226, 277]

# NOW >>> lwed_ind_list (24)
# [280, 269, 242, 236, 251, 265, 261, 267, 243, 230, 253, 241, 252, 234, 260, 291, 256, 258, 245, 273, 228, 227, 231, 277]


simp_ind_list = []
for snp in simp_ind_1:
    if str(snp) in snps_that_arent_working:
        simp_ind_list.append(snp)

#>>> simp_ind_list
#[349, 353, 338, 352, 307, 330, 303, 322, 292, 341, 325, 310, 317, 334, 312, 308, 296, 345, 300, 295, 329, 316, \
#    343, 315, 323, 344, 293, 335, 331, 347, 318, 311, 297, 356]

# NOW >>> simp_ind_list (12)
# [349, 353, 338, 330, 303, 292, 310, 300, 295, 343, 344, 347]

for snp in sex_1:
    if str(snp) in snps_that_arent_working:
        print(snp)

# sex snps = 216, 201, 220

for snp in species_1:
    if str(snp) in snps_that_arent_working:
        print(snp)

# species snps = none

snp_stats="/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/snp_stats.tsv"

# Let's now look at the numbers
import statistics
dict_cov_stats_per_snp = {}
with open(snp_stats, 'w') as output:
    output.write("Set" + "\t" + "snp" + "\t" + "mean_coverage" + "\t" + "mode_coverage" + "\t" + "range_coverage" + "\t" + "stdev_coverage" + "\t" + "variance_coverage" + "\t" + "1st_quantile" + "\t" + "2nd_quantile" + "\t" + "3rd_quantile" + "\n")
    for snp in dict_cov_per_snp:
        mean_cov = statistics.mean(dict_cov_per_snp[snp])
        mode_cov = statistics.mode(dict_cov_per_snp[snp])
        range_cov = max(dict_cov_per_snp[snp]) - min(dict_cov_per_snp[snp])
        stdev_cov = statistics.stdev(dict_cov_per_snp[snp])
        var_cov = statistics.variance(dict_cov_per_snp[snp])
        quantiles_cov = statistics.quantiles(dict_cov_per_snp[snp])
        if snp in ind_id_dir_cov:
            set = "Individual_id_dual"
        if snp in ind_id_lwed_cov:
            set = "LWED_ind_id"
        if snp in ind_id_simp_cov:
            set = "SIMP_ind_id"
        if snp in sex_id_cov:
            set = "Sex_id"
        if snp in species_id_cov:
            set = "Species_id"
        dict_cov_stats_per_snp[snp] = set, mean_cov, mode_cov, range_cov, stdev_cov, var_cov, quantiles_cov[0], quantiles_cov[1], quantiles_cov[2]
        output.write(set + "\t" + str(snp) + "\t" + str(mean_cov) + "\t" + str(mode_cov) + "\t" + str(range_cov) + "\t" + str(stdev_cov) + "\t" + str(var_cov) + "\t" + str(quantiles_cov[0]) + "\t" + str(quantiles_cov[1]) + "\t" + str(quantiles_cov[2]) + "\n")

dict_cov_per_snp
ind_list
coverage_tsv= "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/coverage.tsv"
with open(coverage_tsv, 'w') as output:
    output.write("set" + "\t" + "snp" + "\t" + '\t'.join(ind_list) + "\n")
    for snp in dict_cov_per_snp:
        if snp in ind_id_dir_cov:
            set = "Individual_id_dual"
        if snp in ind_id_lwed_cov:
            set = "LWED_ind_id"
        if snp in ind_id_simp_cov:
            set = "SIMP_ind_id"
        if snp in sex_id_cov:
            set = "Sex_id"
        if snp in species_id_cov:
            set = "Species_id"
        list = dict_cov_per_snp[snp]
        output.write(set + "\t" + str(snp) + "\t" + '\t'.join(str(e) for e in list) + "\n")

# Okay after talking to Katja I know what she wanted me to get: add the sum of coverage for each snp across all individuals after removing the individuals 
# that failed, and let's plot a histogram and look at the distribution.


sum_dict = {}
failed_inds = ["tamOpt-025_S16", "tamOpt-004_S4", "tamOpt-022_S25", "tamOpt-007_S7", "tamOpt-035_S33"]
less_than_280 = []
zero = []
valid = []

for snp in snp_dir:
    snp_sum = 0
    for ind in snp_dir[snp]:
        if ind not in failed_inds:
            sum = int(snp_dir[snp][ind][1]) + int(snp_dir[snp][ind][3])
            snp_sum = snp_sum + sum
    if snp_sum != 0:
        sum_dict[snp] = snp_sum
    if snp_sum == 0:
        zero.append(snp)
    if snp_sum < 280 and snp_sum > 0:
        less_than_280.append(snp)
    if snp_sum != 0 and snp_sum >= 280:
        valid.append(snp)

sum_dict_sorted = sorted(sum_dict.items(), key = lambda x:x[1])
snp_name = []
sum = []
for i in sum_dict_sorted:
    snp_name.append(str(i[0]))
    sum.append(int(i[1]))

plt.bar(snp_name, sum)
plt.xticks(range(0,len(snp_name)),snp_name, rotation = 45, fontsize=6)
plt.show()

plt.bar(snp_name[:146], sum[:146])
plt.xticks(range(0,len(snp_name[:146])),snp_name[:146], rotation = 45, fontsize=6)
plt.axhline(y = 280, color = 'r', linestyle = '-')
plt.show()

plt.bar(snp_name[146:], sum[146:])
plt.xticks(range(0,len(snp_name[146:])),snp_name[146:], rotation = 45, fontsize=6)
plt.show()

# the 34 that are below 280:
plt.bar(snp_name[:38], sum[:38])
plt.xticks(range(0,len(snp_name[:38])),snp_name[:38], rotation = 45, fontsize=6)
plt.show()

# another way to plot these:
less_280_dict = {}
for snp in dict_cov_per_snp:
    if snp in less_than_280:
        less_280_dict[snp] = dict_cov_per_snp[snp]

less_280 = pd.DataFrame.from_dict(less_280_dict)
less_280_reindex = less_280.reindex(less_280.mean().sort_values().index,axis = 1)
list_snps = []
for key in less_280_reindex:
    list_snps.append(str(key))


plt.boxplot(less_280_reindex)
plt.xticks(range(0,len(list_snps)),list_snps, rotation = 45, fontsize=10)
plt.show()

# what about the failed ones!? let's see:
zero_dual = []
zero_lwed = []
zero_simp = []
zero_spp = []
zero_sex = []
for snp in zero:
    if snp in ind_id_dual:
        zero_dual.append(snp)
    elif snp in ind_id_lwed:
        zero_lwed.append(snp)
    elif snp in ind_id_simp:
        zero_simp.append(snp)
    elif snp in sex_id:
        zero_sex.append(snp)
    elif snp in species_id:
        zero_spp.append(snp)



# let's find out what set the less than 280 belong to:
ind_id_dual
ind_id_lwed
ind_id_simp
sex_id
species_id

less_than_280_dual_ind = []
less_than_280_lwed = []
less_than_280_simp = []
less_than_280_sex = []
less_than_280_spp = []

for i in less_than_280:
    if i in ind_id_dual:
        less_than_280_dual_ind.append(i)
    elif i in ind_id_lwed:
        less_than_280_lwed.append(i)
    elif i in ind_id_simp:
        less_than_280_simp.append(i)
    elif i in sex_id:
        less_than_280_sex.append(i)
    elif i in species_id:
        less_than_280_spp.append(i)
    else: 
        print(i)

# let's look at the valid ones 
valid_dual = []
valid_lwed = []
valid_simp = []
valid_spp = []
valid_sex = []

for snp in valid:
    if snp in ind_id_dual:
        valid_dual.append(snp)
    elif snp in ind_id_lwed:
        valid_lwed.append(snp)
    elif snp in ind_id_simp:
        valid_simp.append(snp)
    elif snp in sex_id:
        valid_sex.append(snp)
    elif snp in species_id:
        valid_spp.append(snp)

# now the bp3 error snps:
bp3_error_dual = []
bp3_error_lwed = []
bp3_error_simp = []
bp3_error_spp = []
bp3_error_sex = []

for snp in snps_that_arent_working:
    if int(snp) in ind_id_dual:
        bp3_error_dual.append(snp)
    elif int(snp) in ind_id_lwed:
        bp3_error_lwed.append(snp)
    elif int(snp) in ind_id_simp:
        bp3_error_simp.append(snp)
    elif int(snp) in sex_id:
        bp3_error_sex.append(snp)
    elif int(snp) in species_id:
        bp3_error_spp.append(snp)


# snps that passed bp3 error
pass_bp3_dual = []
pass_bp3_lwed = []
pass_bp3_simp = []
pass_bp3_spp = []
pass_bp3_sex = []


for snp in dict_cov_per_snp:
    if snp in ind_id_dual:
        pass_bp3_dual.append(snp)
    elif snp in ind_id_lwed:
        pass_bp3_lwed.append(snp)
    elif snp in ind_id_simp:
        pass_bp3_simp.append(snp)
    elif snp in sex_id:
        pass_bp3_sex.append(snp)
    elif snp in species_id:
        pass_bp3_spp.append(snp)

# snps that are bp3 errors and are producing snps (MORE THAN 280)

bp3_error_prod_reads_dual = []
bp3_error_prod_reads_lwed = []
bp3_error_prod_reads_simp = []
bp3_error_prod_reads_spp = []
bp3_error_prod_reads_sex = []

for snp in valid_dual:
    if str(snp) in bp3_error_dual:
        bp3_error_prod_reads_dual.append(snp)

for snp in valid_lwed:
    if str(snp) in bp3_error_lwed:
        bp3_error_prod_reads_lwed.append(snp)

for snp in valid_simp:
    if str(snp) in bp3_error_simp:
        bp3_error_prod_reads_simp.append(snp)

for snp in valid_spp:
    if str(snp) in bp3_error_spp:
        bp3_error_prod_reads_spp.append(snp)

for snp in valid_sex:
    if str(snp) in bp3_error_sex:
        bp3_error_prod_reads_sex.append(snp)


no_bp3_error_prod_reads_dual = []
no_bp3_error_prod_reads_lwed = []
no_bp3_error_prod_reads_simp = []
no_bp3_error_prod_reads_spp = []
no_bp3_error_prod_reads_sex = []

for snp in valid_dual:
    if str(snp) not in bp3_error_dual:
        no_bp3_error_prod_reads_dual.append(snp)

for snp in valid_lwed:
    if str(snp) not in bp3_error_lwed:
        no_bp3_error_prod_reads_lwed.append(snp)

for snp in valid_simp:
    if str(snp) not in bp3_error_simp:
        no_bp3_error_prod_reads_simp.append(snp)

for snp in valid_spp:
    if str(snp) not in bp3_error_spp:
        no_bp3_error_prod_reads_spp.append(snp)

for snp in valid_sex:
    if str(snp) not in bp3_error_sex:
        no_bp3_error_prod_reads_sex.append(snp)


# coverage calculation per individual
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt


coverage_tsv= "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/coverage.tsv"
cov_df = pd.read_csv(coverage_tsv, sep = "\t")
new_cov_df=pd.melt(cov_df.drop(["set"],axis = 1).drop(["snp"], axis = 1))
ax = sns.boxplot(x="variable", y="value", data=new_cov_df)
plt.xticks(rotation=45)
ax.figure.tight_layout()



dropped_df= cov_df.drop(["set"],axis = 1).drop(["snp"], axis = 1)
dropped_df.columns = ["73", "644", "326", "466", "NA", "10", "27", "10", "516", "216", "203", "1", "10", "73", "383", "1223", "610", "216", "27", "182", "388", "429", "182", "203", "159", "540", "223", "524", "75", "27", "148", "182", "223"]


coverage_tsv= "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/coverage_plotting_per_sample.tsv"
cov_df = pd.read_csv(coverage_tsv, sep = "\t", header=[0,1])
cov_df_stack=cov_df.stack(level=(0,1)).reset_index(name='value').drop(['level_0'], axis = 1)

# WORKS GREAT
sns.set(rc = {"figure.figsize":(17,10)})
plot= sns.catplot('value', hue='level_2', y='level_1', data=cov_df_stack, kind='box', dodge=False, order=["tamOpt-026_S18", "tamOpt-004_S4", "tamOpt-032_S30", "tamOpt-029_S24", "tamOpt-001_S1", "tamOpt-034_S32", "tamOpt-028_S22", "tamOpt-006_S6", "tamOpt-030_S26", "tamOpt-002_S2", "tamOpt-005_S5", "tamOpt-033_S31", "tamOpt-027_S20", "tamOpt-031_S28", "tamOpt-003_S3", "tamOpt-016_S14", "tamOpt-008_S8", "tamOpt-018_S17", "tamOpt-025_S16", "tamOpt-017_S15", "tamOpt-020_S21", "tamOpt-010_S9", "tamOpt-023_S27","tamOpt-015_S13", "tamOpt-012_S10", "tamOpt-019_S19", "tamOpt-024_S29","tamOpt-022_S25","tamOpt-007_S7","tamOpt-021_S23","tamOpt-013_S11", "tamOpt-014_S12", "tamOpt-035_S33"])
#plt.xlim(0,100)

list=[]
for index in cov_df:
    print(index[0])

# trying some better plots
coverage_tsv= "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/coverage.tsv"
cov_df = pd.read_csv(coverage_tsv, sep = "\t")
new_cov_df = cov_df.rename(columns = {'tamOpt-001_S1': "1",'tamOpt-016_S14': "16", 'tamOpt-008_S8': "8", 'tamOpt-018_S17': "18", 'tamOpt-025_S16':"25", 'tamOpt-032_S30':"32", 'tamOpt-027_S20':27, 'tamOpt-026_S18': "26", \
    'tamOpt-017_S15': "17", 'tamOpt-031_S28': "31", 'tamOpt-002_S2': "2", 'tamOpt-020_S21': "20", 'tamOpt-004_S4': "4", 'tamOpt-029_S24': "29", 'tamOpt-010_S9': "10", 'tamOpt-023_S27': "23", \
        'tamOpt-015_S13': "15", 'tamOpt-003_S3':"3", 'tamOpt-033_S31':"33", 'tamOpt-006_S6':"6",'tamOpt-012_S10':"12", 'tamOpt-019_S19':"19", 'tamOpt-034_S32':"34", 'tamOpt-030_S26':"30", \
            'tamOpt-024_S29':"24", 'tamOpt-022_S25':"22", 'tamOpt-007_S7':"7", 'tamOpt-021_S23':"21",'tamOpt-013_S11':"13", 'tamOpt-005_S5':"5", 'tamOpt-014_S12':"14", 'tamOpt-028_S22':"28", 'tamOpt-035_S33': "35"})
new_cov_df_melt=pd.melt(new_cov_df.drop(["set"],axis = 1).drop(["snp"], axis = 1))
#new_cov_df_melt=pd.melt(new_cov_df)
sns.set(rc = {"figure.figsize":(14,3)})
sns.set(rc = {"figure.figsize":(14,10)})

ax= sns.boxplot(x="variable", y="value", data=new_cov_df_melt, color = "b")
#plt.ylim(-10,300)
#ax.axhline(10, color = 'red')
plt.ylabel("Reads per SNP")
plt.xlabel("Samples")
plt.title("Boxplot of reads obtained per SNP for every sample")
ax.figure.tight_layout()

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/cov_small.png", dpi=300)


sns.set(rc = {"figure.figsize":(14,10)})
ax= sns.boxplot(x="variable", y="value", data=new_cov_df_melt, color = "b")
plt.ylim(-10,300)
ax.axhline(10, color = 'red')
plt.ylabel("Reads per SNP")
plt.xlabel("Samples")
plt.title("Boxplot of reads obtained per SNP for every sample - Zoomed in axis")
ax.figure.tight_layout()

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/cov_zoom.png", dpi=300)

