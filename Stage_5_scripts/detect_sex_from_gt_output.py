# detect sex
import csv
import glob
from traceback import format_exception_only
import pandas as pd
from pandas import DataFrame
from collections import defaultdict

print(glob.glob("/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/updated_genos/*.genos"))
dict_of_inds = {}
for file in glob.glob("/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/updated_genos/*.genos"):
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

sex_id = [195,196,197,198,200,201,203,205,208,209,211,213,215,216,218,219,220,222]

sex_counts = "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/sex_counts_2.tsv"

snp_gen = {'203': {'male': ('G', 'A'), 'female': ('G', 'G')}, '205': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
    '208': {'male': ('G', 'A'), 'female': ('G', 'G')}, '218': {'male': ('C', 'G'), 'female': ('C', 'C')}}

sex_dict = {}

with open(sex_counts, 'w') as output:
    output.write("Sample" + "\t" + "snp" + "\t" + "genotype" + "\t" + "a1_counts" + "\t" + "a2_counts" + "\t" + "a1a2" + "\n")
    for ind in dict_of_inds:
        sex_dict_per_ind = {}
        for snp in dict_of_inds[ind]:
            if snp in [203, 205, 208, 218]:
                if dict_of_inds[ind][snp][5] == "A1HOM":
                    if dict_of_inds[ind][snp][0] == snp_gen[str(snp)]["female"][0]:
                        #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen[str(snp)]["female"][0] + snp_gen[str(snp)]["female"][0] + "\n")
                    elif dict_of_inds[ind][snp][0] == snp_gen[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen[str(snp)]["male"][1] + snp_gen[str(snp)]["male"][1] + "\n")
                if dict_of_inds[ind][snp][5] == "A2HOM":
                    if dict_of_inds[ind][snp][2] == snp_gen[str(snp)]["female"][0]:
                        #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen[str(snp)]["female"][0] + snp_gen[str(snp)]["female"][0] + "\n")
                    elif dict_of_inds[ind][snp][2] == snp_gen[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen[str(snp)]["male"][1] + snp_gen[str(snp)]["male"][1] + "\n")
                if dict_of_inds[ind][snp][5] == "HET":
                    if ((dict_of_inds[ind][snp][0] == snp_gen[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen[str(snp)]['male'][1])):
                        #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen[str(snp)]["male"][0] +  snp_gen[str(snp)]["male"][1] + "\n")
                    #else:
                        #print("This doesn't meet the condition: " + ind + str(snp) + dict_of_inds[ind][snp])
                if dict_of_inds[ind][snp][5] == "nan":
                    #sex_dict_per_ind[snp] = "nan", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                    output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "nan" + "nan" + "\n")
        #sex_dict[ind] = sex_dict_per_ind

# dictionary with the actual genotypes:
sex_genotypes = {'816360': {'male': ('G', 'A'), 'female': ('G', 'G')}, '900829': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
    '913170': {'male': ('G', 'A'), 'female': ('G', 'G')}, '1230640': {'male': ('G', 'A'), 'female': ('G', 'G')}}

snp_gen = {'203': {'male': ('G', 'A'), 'female': ('G', 'G')}, '205': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
    '208': {'male': ('G', 'A'), 'female': ('G', 'G')}, '218': {'male': ('C', 'G'), 'female': ('C', 'C')}}


snp_gen = {'203': {'male': ('G', 'A'), 'female': ('G', 'G')}, '205': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
    '208': {'male': ('G', 'A'), 'female': ('G', 'G')}, '218': {'male': ('C', 'G'), 'female': ('C', 'C')}}


# ____________ trying to include genotypes for all the snps, be they dual or species specific.

import csv
import glob
from traceback import format_exception_only
import pandas as pd
from pandas import DataFrame
from collections import defaultdict

print(glob.glob("/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/updated_genos/*.genos"))
dict_of_inds = {}
for file in glob.glob("/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/updated_genos/*.genos"):
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

sex_id = [195,196,197,198,200,201,203,205,208,209,211,213,215,216,218,219,220,222]

sex_counts = "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/sex_counts_all.tsv"
sex_counts = "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/sex_counts_3.tsv"

simp_genotypes = {'364786': {'male': ('G', 'C'), 'female': ('G', 'G')}, \
    '383098': {'male': ('T', 'C'), 'female': ('C', 'C')}, \
        '790719': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
            '816360': {'male': ('G', 'A'), 'female': ('A', 'A')}, \
                '900829': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                    '916658': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                        '1076256': {'male': ('T', 'A'), 'female': ('T', 'T')}, \
                            '1100033': {'male': ('T', 'C'), 'female': ('T', 'T')}, \
                                '1230640': {'male': ('G', 'A'), 'female': ('A', 'A')}, \
                                    '4369612': {'male': ('G', 'T'), 'female': ('G', 'G')}}

snp_gen_simp = {'196': {'male': ('G', 'C'), 'female': ('G', 'G')}, \
    '198': {'male': ('T', 'C'), 'female': ('C', 'C')}, \
        '200': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
            '203': {'male': ('G', 'A'), 'female': ('A', 'A')}, \
                '205': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                    '209': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                        '213': {'male': ('T', 'A'), 'female': ('T', 'T')}, \
                            '215': {'male': ('T', 'C'), 'female': ('T', 'T')}, \
                                '218': {'male': ('G', 'A'), 'female': ('A', 'A')}, \
                                    '222': {'male': ('G', 'T'), 'female': ('G', 'G')}, \
                                        "197": {'male': ('A', 'G'), 'female': ('A', 'A')}, \
                                            "219": {'male': ('C', 'G'), 'female': ('C', 'C')}} # 197 and 219 have less than 280 reads!


lwed_genotypes={'332842': {'male': ('C', 'T'), 'female': ('T', 'T')}, \
    '364786': {'male': ('G', 'C'), 'female': ('G', 'G')}, \
        '790719': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
            '900829': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                '913170': {'male': ('G', 'A'), 'female': ('G', 'G')}, \
                    '916658': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                        '1076256': {'male': ('T', 'A'), 'female': ('T', 'T')}, \
                            '1100033': {'male': ('T', 'C'), 'female': ('T', 'T')}, \
                                '4369612': {'male': ('G', 'T'), 'female': ('G', 'G')}}

snp_gen_lwed={'195': {'male': ('C', 'T'), 'female': ('T', 'T')}, \
    '196': {'male': ('G', 'C'), 'female': ('G', 'G')}, \
        '200': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
            '205': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                '208': {'male': ('G', 'A'), 'female': ('G', 'G')}, \
                    '209': {'male': ('C', 'T'), 'female': ('C', 'C')}, \
                        '213': {'male': ('T', 'A'), 'female': ('T', 'T')}, \
                            '215': {'male': ('T', 'C'), 'female': ('T', 'T')}, \
                                '222': {'male': ('G', 'T'), 'female': ('G', 'G')}, \
                                    "211": {'male': ('G', 'A'), 'female': ('A', 'A')},\
                                        "219": {'male': ('C', 'G'), 'female': ('C', 'C')}} # 211 and 219 have less than 280 reads!


# the following dictionary will do for now but you need to write a way to get this from the seq data.

ind_spp= {'tamOpt-001_S1': "SIMP", "tamOpt-016_S14": "SIMP", "tamOpt-008_S8": "LWED", "tamOpt-018_S17": "SIMP", "tamOpt-025_S16": "NA",\
            "tamOpt-032_S30": "LWED", "tamOpt-027_S20": "LWED", "tamOpt-026_S18": "LWED", "tamOpt-017_S15": "SIMP","tamOpt-031_S28": "LWED", \
                "tamOpt-002_S2": "LWED", "tamOpt-020_S21": "SIMP", "tamOpt-004_S4": "NA", "tamOpt-029_S24": "SIMP", "tamOpt-010_S9": "LWED",\
                    "tamOpt-023_S27": "SIMP", "tamOpt-015_S13": "SIMP", "tamOpt-003_S3": "LWED", "tamOpt-033_S31": "LWED", "tamOpt-006_S6": "SIMP", \
                        "tamOpt-012_S10": "LWED", "tamOpt-019_S19": "SIMP", "tamOpt-034_S32": "SIMP", "tamOpt-030_S26": "LWED", "tamOpt-024_S29": "SIMP",\
                            "tamOpt-022_S25": "NA", "tamOpt-007_S7": "NA", "tamOpt-021_S23": "SIMP", "tamOpt-013_S11": "LWED", "tamOpt-005_S5": "LWED",\
                                "tamOpt-014_S12": "LWED","tamOpt-028_S22": "SIMP","tamOpt-035_S33": "LWED"}

sex_dict = {}

with open(sex_counts, 'w') as output:
    output.write("Sample" + "\t" + "snp" + "\t" + "genotype" + "\t" + "a1_counts" + "\t" + "a2_counts" + "\t" + "a1a2" + "\t" + "inferred_spp" + "\n")
    for ind in dict_of_inds:
        sex_dict_per_ind = {}
        for snp in dict_of_inds[ind]:
            spp = ind_spp[ind]
            if snp in [208, 195, 211]:
                if spp == "LWED":
                #if snp in [208, 211, 205, 219]:    
                    if dict_of_inds[ind][snp][5] == "A1HOM":
                        if dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["female"][0] + snp_gen_lwed[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["male"][1] + snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "A2HOM":
                        if dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_lwed[str(snp)]["female"][0] + snp_gen_lwed[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["male"][1] + snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "HET":
                        if ((dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]['male'][1])):
                            #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_lwed[str(snp)]["male"][0] +  snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                        #else:
                            #print("This doesn't meet the condition: " + ind + str(snp) + dict_of_inds[ind][snp])
                    elif dict_of_inds[ind][snp][5] == "nan":
                        output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")
                elif spp == "SIMP" and (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3]) > 0) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
                elif spp == "NA" and (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3]) > 0) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_NA" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
            elif snp in [198, 203, 218, 197]:
                if spp == "SIMP":
                    if dict_of_inds[ind][snp][5] == "A1HOM":
                        if dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] +"\t" + spp +  "\n")
                        elif dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "A2HOM":
                        if dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "HET":
                        if ((dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][1])):
                            #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["male"][0] +  snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "nan":
                        output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")
                elif spp == "LWED" and (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3]) > 0) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
                elif spp == "NA" and (int(dict_of_inds[ind][snp][1]) + int(dict_of_inds[ind][snp][3]) >0) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_NA" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp +  "\n")
            elif snp in [200, 205, 209, 215, 222, 219]:
                if dict_of_inds[ind][snp][5] == "A1HOM":
                    if dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["female"][0]:
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp +  "\n")
                    elif dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "A2HOM":
                    if dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["female"][0]:
                        #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "HET":
                    if ((dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][1])):
                        #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["male"][0] +  snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "nan":
                    output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")


# We also have these variants: [197, 211, 219] which are producing less than 280 reads across all individuals but that we might as well use since they are for sex
# let's first check that they don't have the bp3 error. (they don't!)

# 197 is informative in SIMPs only
# 211 is informative in LWEDs only
# 219 is informative in both

# included!


# trying to be more stringent - just exploring here

sex_counts = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/sex_counts_het10s.tsv"


with open(sex_counts, 'w') as output:
    output.write("Sample" + "\t" + "snp" + "\t" + "genotype" + "\t" + "a1_counts" + "\t" + "a2_counts" + "\t" + "a1a2" + "\t" + "inferred_spp" + "\n")
    for ind in dict_of_inds:
        sex_dict_per_ind = {}
        for snp in dict_of_inds[ind]:
            spp = ind_spp[ind]
            if snp in [208, 195, 211]:
                if spp == "LWED" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10):  
                    if dict_of_inds[ind][snp][5] == "A1HOM":
                        if dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["female"][0] + snp_gen_lwed[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["male"][1] + snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "A2HOM":
                        if dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_lwed[str(snp)]["female"][0] + snp_gen_lwed[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_lwed[str(snp)]["male"][1] + snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "HET":
                        if (int(dict_of_inds[ind][snp][1]) >= 10) and (int(dict_of_inds[ind][snp][3]) >=10):
                            if ((dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_lwed[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_lwed[str(snp)]['male'][1])):
                                #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                                output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_lwed[str(snp)]["male"][0] +  snp_gen_lwed[str(snp)]["male"][1] + "\t" + spp + "\n")
                            #else:
                                #print("This doesn't meet the condition: " + ind + str(snp) + dict_of_inds[ind][snp])
                    elif dict_of_inds[ind][snp][5] == "nan":
                        output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")
                elif spp == "SIMP" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
                elif spp == "NA" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10):
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_NA" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
            elif snp in [198, 203, 218, 197]:
                if spp == "SIMP" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10):
                    if dict_of_inds[ind][snp][5] == "A1HOM":
                        if dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] +"\t" + spp +  "\n")
                        elif dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "A2HOM":
                        if dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["female"][0]:
                            #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp + "\n")
                        elif dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["male"][1]:
                            output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "HET":
                        if (int(dict_of_inds[ind][snp][1]) >= 10) and (int(dict_of_inds[ind][snp][3]) >= 10):
                            if ((dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][1])):
                                #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                                output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous_simp" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["male"][0] +  snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][5] == "nan":
                        output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")
                elif spp == "LWED" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10) :
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_lwed" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp + "\n")
                elif spp == "NA" and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10):
                    output.write(ind + "\t" + str(snp) + "\t" + "reads_present_but_not_informative_for_NA" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + str(dict_of_inds[ind][snp][0]) + str(dict_of_inds[ind][snp][2]) + "\t" + spp +  "\n")
            elif snp in [200, 205, 209, 215, 222, 219] and (int(dict_of_inds[ind][snp][1]) >= 10 or int(dict_of_inds[ind][snp][3]) >= 10):
                if dict_of_inds[ind][snp][5] == "A1HOM":
                    if dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["female"][0]:
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp +  "\n")
                    elif dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "A2HOM":
                    if dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["female"][0]:
                        #sex_dict_per_ind[snp] = "Female", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "X_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["female"][0] + snp_gen_simp[str(snp)]["female"][0] + "\t" + spp + "\n")
                    elif dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]["male"][1]:
                        output.write(ind + "\t" + str(snp) + "\t" + "Y_chrom_homozyg" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" + snp_gen_simp[str(snp)]["male"][1] + snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "HET":
                    if (int(dict_of_inds[ind][snp][1]) >= 10) and (int(dict_of_inds[ind][snp][3]) >=10):
                        if ((dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][1])) or ((dict_of_inds[ind][snp][2] == snp_gen_simp[str(snp)]['male'][0]) and (dict_of_inds[ind][snp][0] == snp_gen_simp[str(snp)]['male'][1])):
                            #sex_dict_per_ind[snp] = "Male", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                            output.write(ind + "\t" + str(snp) + "\t" + "Heterozygous" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\t" +  snp_gen_simp[str(snp)]["male"][0] +  snp_gen_simp[str(snp)]["male"][1] + "\t" + spp + "\n")
                elif dict_of_inds[ind][snp][5] == "nan":
                    output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\n")


