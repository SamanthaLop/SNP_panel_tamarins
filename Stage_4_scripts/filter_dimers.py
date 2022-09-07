# remember to load bioinfo-tools python/3.8.7 and biopython/1.78

# import modules
import numpy as np 
import pandas as pd
import argparse

#dimer_file="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/ANGSD_consensus_seqs/split_mfeoutput/dedup.101462_127.bam_angsd_CM018917_dimerlist_3.tsv"
#dimer_list=pd.read_csv(dimer_file, sep='\t')

# define the function
def filter_df_tm(dimer_list, output_tsv, lower_limit_tm, upper_limit_tm):
    dimer_list=pd.read_csv(dimer_list, sep='\t')
    filtered_df=dimer_list[(dimer_list.Tm_C > lower_limit_tm) & (dimer_list.Tm_C < upper_limit_tm)]
    filtered_df.to_csv(output_tsv,sep="\t")

# parse arguments
parser = argparse.ArgumentParser(description="Function to filter dimer list tsv from MFEprimer output")
parser.add_argument('-i', '--input', type=str, help='Input tsv dimer list', required=True)
parser.add_argument('-o', '--output', type=str, help='Output filtered tsv dimer list', required=True)
parser.add_argument('-u', '--uppertm', type=int, help='Upper tm limit integer in C, defaults to 65')
parser.add_argument('-l', '--lowertm', type=int, help ='Lower tm limit integer in C, defaults to 55')
args= parser.parse_args()

dimer_list=args.input
output_tsv=args.output

if args.uppertm:
    upper_limit_tm=args.uppertm
else:
    upper_limit_tm=65
if args.lowertm:
    lower_limit_tm=args.lowertm
else:
    lower_limit_tm=55

filter_df_tm(dimer_list, output_tsv, lower_limit_tm,upper_limit_tm)
