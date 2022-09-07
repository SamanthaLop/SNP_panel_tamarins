from pysam import VariantFile
import argparse
import os
import sys
import csv

# This script will output a tsv file with four columns, <SNP> <start> <stop> <length>. It looks for consecutive regions above a certain input threshold within the input vcf
# It will also look in the consecutive region for a SNP from our input SNP list and output the SNP, the start and stop of the consecutive region, and the length of the consecutive
# region. The input list of SNPs is assumed to be composed of two columns separated by tabs: <chrom> <position>

# This script is run like this: python3 consecutive_positions_in_vcf.py -i <input.vcf.gz> -l <list_of_snps.tsv> -t <threshold (int)> -o <output.tsv>

# create a position list from all positions present in our vcf file
pos_list=[]
def parse_vcf(input_vcf):
    with VariantFile(input_vcf) as vcf:
        for rec in vcf.fetch():
            pos_list.append(rec.pos)
    

# define the function to find the ranges of consecutive numbers in the list
def ranges(nums):
    nums= sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums [-1:])
    return list(zip(edges, edges))

# create a snp_list with all the snps to check if our consecutive region contains a snp
snp_list=[]
def parse_snp_list(input_list):
    with open(input_list) as SNPsf:
        rd = csv.reader(SNPsf, delimiter="\t")
        for row in rd:
            snp_list.append(row[1])

# use the input threshold, ranges_list and snp_list to add SNPs and their consecutive site to our output if they meet our requirements
def find_consecutive_sites(output, ranges_list, snp_list, threshold):
    with open(output,"w") as out:
        out.write(str("SNP" + "\t" + str("start") + "\t" + str("stop") + "\t" + str("length" + "\n")))
        for i in ranges_list:
            diff=i[1]-i[0]
            if diff > threshold:
                for SNP in snp_list:
                    if int(SNP) in range(i[0], i[1]):
                        out.write(str(SNP) + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(diff) + "\n")

# define our input arguments
parser = argparse.ArgumentParser(description="Find consecutive sites in a filtered VCF")

parser.add_argument('-i', '--input-vcf', type=str, help='Input vcf (can be gzipped)', required=True)
parser.add_argument('-l', '--input-list', type=str, help='Input list of SNPs in tsv (chr \t pos)', required=True)
parser.add_argument('-t', '--threshold', type=int, help="Input minimum number of consecutive sites", required=True)
parser.add_argument('-o', '--output', type=str, help='Output tsv file.')

args = parser.parse_args()

input_vcf = args.input_vcf
input_list = args.input_list
threshold = args.threshold
output = args.output

# run our functions
parse_vcf(input_vcf)

ranges_list = ranges(pos_list)

parse_snp_list(input_list)

output = find_consecutive_sites(output, ranges_list, snp_list, threshold)




#__________________initial code written, without adapting to run with functions and input, saving for my records:
#input_vcf="/proj/proj_name/nobackup/SAM/flanking_regions/female_het_SNPs/700_flanking_intervals_vcfs/gatk_filtered_female_het_700_incl_SNP_chromname_700_flanking.vcf.gz.vcf_repet$
# make the list of positions:
#with VariantFile(input_vcf) as vcf:
#    pos_list=[]
#    for rec in vcf.fetch():
#        pos_list.append(rec.pos)

# feed our data (list of positions) to the function we just defined
#ranges_list=ranges(pos_list)

# make a list of SNPs
#SNPs_file="/proj/proj_name/nobackup/SAM/flanking_regions/female_het_SNPs/lists/female_het_SNP_list.tsv"
#with open(input_list) as SNPsf:
#    snp_list=[]
#    rd = csv.reader(SNPsf, delimiter="\t")
#    for row in rd:
#        snp_list.append(row[1])

# look for the flanking regions that span 350 positions that include our SNP
#output="/proj/proj_name/nobackup/SAM/flanking_regions/female_het_SNPs/700_flanking_intervals_vcfs/consecutive_regions.tsv"
#with open(output,"w") as out:
#    out.write(str("SNP" + "\t" + str("start") + "\t" + str("stop") + "\t" + str("length" + "\n")))
#    for i in ranges_list:
#        diff=i[1]-i[0]
#        if diff > 350:
#            for SNP in snp_list:
#                if int(SNP) in range(i[0], i[1]):
#                    out.write(str(SNP) + "\t" + str(i[0]) + "\t" + str(i[1]) + "\t" + str(diff) + "\n")
