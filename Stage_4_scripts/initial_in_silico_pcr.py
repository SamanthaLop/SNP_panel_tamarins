# import modules
import numpy as np
import pandas as pd
import argparse

# This script contains the exploration of the results of the first in silico PCR that we ran, it contains comments with the results of certain commands that were run. 
# It is meant as an exploratory script to be run one command at a time on an interactive python terminal.

# input file
amplicon_details_file="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/ANGSD_consensus_seqs/REAL_TEST/101479/MFE_output_101479_real_primers_amplicon_merged_4_5.tsv"

# parse input into dataframe
amplicon_data_frame=pd.read_csv(amplicon_details_file, sep='\t')
#calculate number of amplicons
initial_number_of_amplicons=len(amplicon_data_frame.Amp_ID) # 528926

# filter out <53 tm and >64 tm (both forward and reverse)
# this tm was calculated by looking at the max MFEprimer-calculated Tm (53.87) and the min MFEprimer calculated Tm (63.18) (obtained from the primerlist).
ampl_Tm55_67=amplicon_data_frame[(amplicon_data_frame.FpTm_C >53) & (amplicon_data_frame.FpTm_C <64) & (amplicon_data_frame.RpTm_C >53) & (amplicon_data_frame.RpTm_C <64) ]
# calculate number of amplicons left
number_amplicons_tmfilt=len(ampl_Tm55_67.Amp_ID) # 1961

# filter out >1000bp length
ampl_Tm55_67_Lenless1000=ampl_Tm55_67[ampl_Tm55_67.Size_bp < 1000]
# calculate number of amplicons left
number_amplicons_sizefilt=len(ampl_Tm55_67_Lenless1000.Amp_ID) # 812

# look at the [:#] number, that's when you stop seeing numbers < 2
top_primers_all_p1=ampl_Tm55_67_Lenless1000["P_1"].value_counts()[:9]
top_primers_all_list_p1=ampl_Tm55_67_Lenless1000["P_1"].value_counts()[:9].index.tolist()
top_primers_all_p2=ampl_Tm55_67_Lenless1000["P_2"].value_counts()[:13]
top_primers_all_list_p2=ampl_Tm55_67_Lenless1000["P_2"].value_counts()[:13].index.tolist()
print("These are the most repeated primers in our full dataframe:\n", "P_1", top_primers_all_list_p1, "P_2", top_primers_all_list_p2)

# These are the most repeated primers in our full dataframe:
# P_1 [364, 363, 380, 379, 120, 359, 119, 553, 99] P_2 [364, 379, 120, 363, 582, 360, 62, 532, 100, 468, 380, 8, 546]

# now let's merge these two lists:

common_primer_list = list(top_primers_all_list_p1)
common_primer_list.extend(x for x in top_primers_all_list_p2 if x not in common_primer_list)
len(common_primer_list) # 17

# of these 17 common primers, let's see what SNPs those are for.
non_specific_species=[]
non_specific_ind_id=[]
non_specific_sex=[]
non_specific_lwed=[]
non_specific_simp=[]
for i in common_primer_list:
    if i < 41:
        non_specific_species.append(i)
    elif i > 40 & i < 345:
        non_specific_ind_id.append(i)
    elif i > 344 & i < 395:
        non_specific_sex.append(i)
    elif i > 394 & i < 527:
        non_specific_lwed.append(i)
    elif i > 526:
        non_specific_simp.append(i)

print("There are", len(non_specific_species), "primers that are not specific (i.e. making many amplicons) for species.\n", \
    "There are", len(non_specific_ind_id), "primers that are not specific (i.e. making many amplicons) for ind_id.\n", \
    "There are", len(non_specific_sex), "primers that are not specific (i.e. making many amplicons) for sex.\n", \
    "There are", len(non_specific_lwed), "primers that are not specific (i.e. making many amplicons) for LWED ind id.\n", \
    "There are", len(non_specific_simp), "primers that are not specific (i.e. making many amplicons) for SIMP ind id.")

# There are 1 primers that are not specific (i.e. making many amplicons) for species.
# There are 16 primers that are not specific (i.e. making many amplicons) for ind_id.
# There are 0 primers that are not specific (i.e. making many amplicons) for sex.
# There are 0 primers that are not specific (i.e. making many amplicons) for LWED ind id.
# There are 0 primers that are not specific (i.e. making many amplicons) for SIMP ind id.

# from the df where we filtered for tm and length, now filter OUT those rows that contain the non-specific primers (the primers that produce more
# than one amplicon)
df_without_non_specific_primers=ampl_Tm55_67_Lenless1000[~(ampl_Tm55_67_Lenless1000.P_1.isin(common_primer_list) | ampl_Tm55_67_Lenless1000.P_2.isin(common_primer_list))]
len(df_without_non_specific_primers) # 214

# NOW we have 241 amplicons left over, 241 plus the 12 SNPs that we should not see due to removing the repetitive primers gives us 253. This means we have
# 76 primers [241 + 12 - 329 = 76] that are not amplifying their targeted regions under our conditions (Tm 53-64)

# create a list of odd numbers, each SNP has a forward (odd number) and reverse (even number). If we can check what odd numbers (forwards) are absent from 
# our amplicon dataframe, we can find out what primers aren't working, and therefore, what primer pairs aren't working.
odd_numbers=[]
for a in range(1, 658, 2):
    odd_numbers.append(a)

# lets see what primers aren't amplifying from the initial df (without filtering the non-specific primers)
valores_P_1=df_without_non_specific_primers["P_1"].values
primers_not_working=[]
for i in odd_numbers:
    if i not in valores_P_1:
        primers_not_working.append(i)

len(primers_not_working) # 91

# NOW 91 primers are not working., minus the 12 repetitive pairs, thats 79. There are three extra ones (91 - 12 - 76 = 3). 

targets=df_without_non_specific_primers[((df_without_non_specific_primers.P_2 % 2) == 0) & (df_without_non_specific_primers.P_1 == df_without_non_specific_primers.P_2 -1)]
number_amplicons_targets=len(targets.Amp_ID) # 238
non_targets=df_without_non_specific_primers[~((df_without_non_specific_primers.P_2 % 2) == 0) & ~(df_without_non_specific_primers.P_1 == df_without_non_specific_primers.P_2 -1)]
number_non_targets=len(non_targets.Amp_ID)  # 3

# NOW we get 238 amplicons that are targets and three that are non targets.

# okay, so now lets see, these primers should be removed given that they amplify a region that is not our target (forward is acting as reverse)
primers_to_be_removed_2=[368,367,84,83,396,395]
df_without_non_specific_primers=ampl_Tm55_67_Lenless1000[~(ampl_Tm55_67_Lenless1000.P_1.isin(common_primer_list) | ampl_Tm55_67_Lenless1000.P_2.isin(common_primer_list))]

# by removing these we effectively remove six amplicons
removed=df_without_non_specific_primers[df_without_non_specific_primers.P_1.isin(primers_to_be_removed_2) | df_without_non_specific_primers.P_2.isin(primers_to_be_removed_2)]
len(removed) # 6, 3 ARE TARGETS AND 3 AREN'T

targets_1=df_without_non_specific_primers[~(df_without_non_specific_primers.P_1.isin(primers_to_be_removed_2) | df_without_non_specific_primers.P_2.isin(primers_to_be_removed_2))]
len(targets_1) # 235

# NOW 235 targets 

# NOW math:
# 329
#- 12  # sites with one or more of the 17 common primers
#____
# 317 
#-  3  # sites F-R # 3 targets, 3 non targets
#____
# 314
#-235  # targets
#____ 
# 79   

# let's try to get these in a list:
common_primer_list_f_and_r=[363, 364, 363, 380, 379, 119, 120, 359, 360, 119, 553, 554, 99, 581, 582, 360, 531, 61, 467, 7, 545, 62, 532, 100, 468, 8, 546]

values_P_1=targets_1["P_1"].values 
passed=[]
not_passed=[]
for i in odd_numbers:
    if i in values_P_1:
        passed.append(i)    # 235
    else:
        not_passed.append(i)   # 94

# NOW length of passed 235, length of not passed = 94. minus the 12 SNPs = 82, minus the RFFR = 79

not_working_primers=[]
for i in not_passed:
    if i not in common_primer_list_f_and_r:
        not_working_primers.append(i)    # now len 82

not_working_primers_1=[]
for i in not_working_primers:
    if i not in primers_to_be_removed_2:
        not_working_primers_1.append(i)    # now len 79

# Finally! we have our final list of 79 forwards (and thus 79 SNPs) that aren't working. Let's see what set these are for.

# but first let's make a list of all the primers that need replacing (the really common ones, the ones not working and the RFFR ones)
# common primers:
common_primer_list_f_and_r
# not working primers + RFFR
not_working_primers
not_working_primers_f_and_r=[]
for i in not_working_primers:
    not_working_primers_f_and_r.append(i)
    not_working_primers_f_and_r.append(i+1)

primers_need_replacing= common_primer_list_f_and_r + not_working_primers_f_and_r
primers_need_replacing.sort()

# let's make sure no replicates

primers_need_replacing = list(dict.fromkeys(primers_need_replacing))

need_replace_species=[]
need_replace_ind_id=[]
need_replace_sex=[]
need_replace_lwed=[]
need_replace_simp=[]
for i in primers_need_replacing:
    if i < 41:
        need_replace_species.append(i)
    elif i > 40 & i < 345:
        need_replace_ind_id.append(i)
    elif i > 344 & i < 395:
        need_replace_sex.append(i)
    elif i > 394 & i < 527:
        need_replace_lwed.append(i)
    elif i > 526:
        need_replace_simp.append(i)

print("There are", len(need_replace_species) / 2, "primer pairs are not working for species.\n", \
    "There are", len(need_replace_ind_id) / 2, "primer pairs that are not working for ind_id.\n", \
    "There are", len(need_replace_sex) / 2, "primers pairs that are not working for sex.\n", \
    "There are", len(need_replace_lwed) / 2, "primers pairs that are not working for LWED ind id.\n", \
    "There are", len(need_replace_simp) / 2, "primers pairs that are not working for SIMP ind id.")

# There are 3.0 primer pairs are not working for species.
# There are 0.0 primers pairs that are not working for sex.
# There are 91.0 primer pairs that are not working for ind_id.
# There are 0.0 primers pairs that are not working for LWED ind id.
# There are 0.0 primers pairs that are not working for SIMP ind id.


# Rachel was kind enough to find replacements for these primers. 
