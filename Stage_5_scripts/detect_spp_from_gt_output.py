# detect spp
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



species_snp_output="/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/chrom_pos_final_snps_final.tsv"
species_counts = "/proj/proj_name/nobackup/SAM/genotyping/third_geno_run/spp_counts.tsv"
spp_ref = pd.read_csv(species_snp_output, sep = "\t")
species_id = [1,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
species_id_12_snps = [1, 5, 8, 9, 12, 13, 14, 15, 16, 17, 19, 20]


spp_dict = {}
with open(species_counts, 'w') as output:
    output.write("Sample" + "\t" + "snp" + "\t" + "genotype" + "\t" + "a1_counts" + "\t" + "a2_counts" + "\n")
    for ind in dict_of_inds:
        spp_dict_per_ind = {}
        for snp in dict_of_inds[ind]:
            if snp in species_id_12_snps:
                for index,row in spp_ref.iterrows():
                    if row["bp3_index"] == snp:
                        simp_var = row["SIMP_variant"]
                        lwed_var = row["LWED_variant"]
                if dict_of_inds[ind][snp][5] == "A1HOM":
                    if (dict_of_inds[ind][snp][0] == simp_var):
                        spp_dict_per_ind[snp] = "SIMP", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "SIMP" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")
                    elif (dict_of_inds[ind][snp][0] == lwed_var):
                        spp_dict_per_ind[snp] = "LWED", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "LWED" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")
                elif dict_of_inds[ind][snp][5] == "A2HOM":
                    if (dict_of_inds[ind][snp][2] == simp_var) :
                        spp_dict_per_ind[snp] = "SIMP", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "SIMP" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")
                    elif (dict_of_inds[ind][snp][2] == lwed_var):
                        spp_dict_per_ind[snp] = "LWED", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                        output.write(ind + "\t" + str(snp) + "\t" + "LWED" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")
                elif dict_of_inds[ind][snp][5] == "HET":
                    spp_dict_per_ind[snp] = "het", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                    output.write(ind + "\t" + str(snp) + "\t" + "het" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")                    
                elif dict_of_inds[ind][snp][5] == "nan":
                    spp_dict_per_ind[snp] = "nan", dict_of_inds[ind][snp][1], dict_of_inds[ind][snp][3]
                    output.write(ind + "\t" + str(snp) + "\t" + "nan" + "\t" + str(dict_of_inds[ind][snp][1]) + "\t" + str(dict_of_inds[ind][snp][3]) + "\n")


