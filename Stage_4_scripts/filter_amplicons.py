# import modules
import numpy as np 
import pandas as pd
import argparse

#python script to filter amplicon description tsv file MFEprimer output with default or chosen Tm and size limits.
#python3 filter_amplicons.py -i <input tsv amplicon description file> -o <output filtered amplicon details file> 
#which will use 65 as the upper Tm limit, 55 as the lower Tm limit and 1000bp as the upper size limit. The script can receive input however, if these need to change:
#python3 filter_amplicons.py -i <input tsv amplicon description file> -o <output filtered amplicon description file> -u <int> -l <int> -s <int>


# define the function
def filter_df_tm(amplicon_details_file, output_tsv, lower_limit_tm, upper_limit_tm, upper_size_limit):
    amplicon_details_file=pd.read_csv(amplicon_details_file, sep='\t')
    filtered_df=amplicon_details_file[(amplicon_details_file.FpTm_C > lower_limit_tm) & (amplicon_details_file.FpTm_C < upper_limit_tm) \
    & (amplicon_details_file.RpTm_C > lower_limit_tm) & (amplicon_details_file.RpTm_C < upper_limit_tm) & (amplicon_details_file.Size_bp < upper_size_limit)]
    filtered_df.to_csv(output_tsv,sep="\t")

# parse arguments
parser = argparse.ArgumentParser(description="Function to filter amplicon details tsv from MFEprimer output")
parser.add_argument('-i', '--input', type=str, help='Input tsv amplicon details', required=True)
parser.add_argument('-o', '--output', type=str, help='Output filtered tsv amplicon details', required=True)
parser.add_argument('-u', '--uppertm', type=int, help='Upper tm limit integer in C, defaults to 65')
parser.add_argument('-l', '--lowertm', type=int, help ='Lower tm limit integer in C, defaults to 55')
parser.add_argument('-s', '--sizeupper', type=int, help='Upper size limit in bp, defaults to 1000')
args= parser.parse_args()

# parse command inputs
amplicon_details_file=args.input
output_tsv=args.output

# optional command inputs, else set to default.
if args.uppertm:
    upper_limit_tm=args.uppertm
else:
    upper_limit_tm=65
if args.lowertm:
    lower_limit_tm=args.lowertm
else:
    lower_limit_tm=55
if args.sizeupper:
    upper_size_limit=args.sizeupper
else:
    upper_size_limit=1000


# run the function
filter_df_tm(amplicon_details_file, output_tsv, lower_limit_tm, upper_limit_tm, upper_size_limit)
