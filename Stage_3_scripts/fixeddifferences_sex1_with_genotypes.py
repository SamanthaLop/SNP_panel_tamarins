from pysam import VariantFile
import argparse
import os
import sys


# this script should be run like this: sbatch fixeddifferences.py -i input.vcf --popfile popfile -o output.tsv


def parse_popfile(popfile):
    population_ids= [] # empty list to hold the pop-ids
    populations = {} # empty dictionary to have a list of what samples belong to which pop
    with open(popfile) as p: 
        for line in p.readlines(): # iterate through the popfile, line by line
            try: # put a try statement here, to allow for samples not to be assigned to a pop. If the statement within this fails, it will simply output that this sample wasn't assigned to a pop and then continue. Without a try/except statement it would have exited the script.
                sample = line.split("\t")[0] # split the line by tab, first item [0] will be the sample
                pop = line.split("\t")[1].rstrip() # second will be the pop assignment, rstrip tells it to remove the trailing newline. Will throw error if a sample has no pop
                if not pop in population_ids: # if the pop name isn't already in the list of ids, we add it
                    if not pop == "": # sometimes the last line without any real text will be added, this just makes sure that no empty lines are added to the pop ids
                        population_ids.append(pop)
                        populations[pop]=[sample] # if this is the first item for the population, add a key with the pop-id to the dictionary, and a list as a value where we add the first sample
                    else: # don't think this does anythin - if pop is empty it shouldnt in principle get here, but in that case it won't assign a sample to an empty pop name
                        print("sample " + sample.rstrip() + " not assigned to a population.")
                else: # if it's not the first time a population appears, append this sample to the value list belonging to the pop-key in the dictionary
                    populations[pop].append(sample)
            except:
                print("sample " + sample.rstrip() + " not assigned to a population")
    return population_ids, populations #(this is for when you're defining a function)

input_vcf="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_23_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz"
popfile="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/popfile"
output = "/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/new_fixed_diff_sex1_genotypes.tsv"

sex_pos = ["332842", "364786", "370770", "383098", "790719", "803345", "816360", "900829", "913170", "916658", "925487", "1076256", "1100033", "1196778", "1230640", "1232589", "1260827", "4369612"]
sex_pos_bp3 = ["195", "196", "197", "198", "200", "201", "203", "205", "208", "209", "211", "213", "215", "216", "218", "219", "220", "222"]

# THIRD geno run, let's get the 12 sex SNPs
#[195, 196, 198, 200, 203, 205, 208, 209, 213, 215, 218, 222]
sex_pos_3 = ["332842", "364786", "383098", "790719", "816360", "900829", "913170", "916658", "1076256", "1100033", "1230640", "4369612"]
sex_pos_4 = ["332842", "364786", "383098", "790719", "816360", "900829", "913170", "916658", "1076256", "1100033", "1230640", "4369612", "370770", "925487", "1232589"] # adding the variants that have <280 reads


population_ids, populations = parse_popfile(popfile)

#fixed_differences = [] # empty list where we will append one dictionary for each position, containing all the relevant information
#def find_fixed_differences(input_vcf,populations,population_ids, maxmissing=None):

fixed_differences = []
genotypes = {}
with VariantFile(input_vcf) as vcf: # parsing the vcf file with pysam
    for snp in sex_pos_4:
        for rec in vcf.fetch(): # fetch() creates an iterator where we can go trhough all the records in the vcf file
            male_is_het = [] # move empty lists here !!!!!
            female_is_hom = []
            female_alleles= []
            female_alleles_1 = []
            male_alleles_1 = []
            sex_dir = {}
            position = rec.pos # will give us the position
            chrom = rec.chrom # the chromosome/scaffold
            rec_alleles = rec.alleles # this will be tuble with all the alleles present at a site, where the ref is the first [0]. For bi-allelic sites, this should always be a length of 2-
            x_fixed= False
            alleles = {pop: [] for pop in population_ids} # a disctionary to hold all the alleles that are in a pop, using dictionary comprehension. It's a way of looping in just one line, could just as well have been created in a loop. 
            if snp == str(position):
                for pop in population_ids:
                    for sample in populations[pop]: # go through each sample in both populations, one at a time
                        a1 = rec.samples[sample]['GT'][0] # will extract first allele that is called in the sample
                        a2 = rec.samples[sample]['GT'][1] # and the second
                        if pop == "female":     # within females
                            if a1 == a2 and not None in [a1,a2]:  # if homozygous (fixed in X)
                                female_is_hom.append(sample)
                                female_alleles.append(a1)
                                female_alleles_1.append(rec_alleles[0])
                                female_alleles_1.append(rec_alleles[1])
                                sex_dir["female"] = rec_alleles[a1], rec_alleles[a2]
                        elif pop == "male":  # within males
                            if not a1 == a2 and not None in [a1,a2]: # if heterozygous (fixed in Y)
                                male_is_het.append(sample)
                                male_alleles_1.append(rec_alleles[0])
                                male_alleles_1.append(rec_alleles[1])
                                sex_dir["male"] = rec_alleles[a1], rec_alleles[a2]
                #genotypes[snp] = sex_dir # I only retrieve the genotypes if I add this step here... Figured it out! look below, replacement SNPs are spp specific
                pop2= population_ids[1]
                pop1= population_ids[0]
                if len(female_is_hom) == len(populations["female"]):
                    if len(set(female_alleles)) == 1:
                        x_fixed = True
                if len(male_is_het) == len(populations["male"]):
                    if x_fixed:
                        fixed_differences.append({"chromosome": chrom, "position": position})
                        genotypes[snp] = sex_dir # only one of the four snps passes this condition.


parser = argparse.ArgumentParser(description="Extract SNPs that are homozygous in females and heterozygous in males")

parser.add_argument('-i', '--input-vcf', type=str, help='Input vcf (can be gzipped)', required=True)
parser.add_argument('--popfile',type=str, help="File with sample in first column and sex - female or male - in the second, separated by a tab")
parser.add_argument('-o', '--output', type=str, help='Output tsv file.')

args = parser.parse_args()

# parse input arguments
input_vcf = args.input_vcf
output = args.output
popfile = args.popfile

# parse popfile
population_ids, populations = parse_popfile(popfile)

# check that the number of sexes is 2:
if not len(population_ids) == 2:
    print("Number of sexes must be two.")
    sys.exit()

# identify sites
fixed_differences=find_fixed_differences(input_vcf, populations, population_ids, maxmissing=None)

# write to output
with open(output,"w") as of:
    of.write("Chromosome\tPosition\n")
    for diff in fixed_differences:
        of.write(str(diff["chromosome"]) + "\t" + str(diff["position"]) + "\n")
    of.close()
