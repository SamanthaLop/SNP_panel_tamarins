from pysam import VariantFile
import argparse
import os
import sys

# Can run it with commandline arguments like python3 fixedDifferences.py -i /your/vcf.vcf -o output.tsv --popfile sample_species.txt
# popfile should be a file listing all the samples in the first column, and which of the two species it belongs to in the second (tab-separated)

def parse_popfile(popfile):
	population_ids = []
	populations = {}
	with open(popfile) as p:
		for line in p.readlines():
			try:
				sample = line.split("\t")[0]
				pop = line.split("\t")[1].rstrip()
				if not pop in population_ids:
					if not pop == "":
						population_ids.append(pop)
						populations[pop] = [sample]
					else:
						print("sample " + sample.rstrip() + " not assigned to a population.")
				else:
					populations[pop].append(sample)
			except:
				print("sample " + sample.rstrip() + " not assigned to a population.")
	return population_ids, populations

def find_fixed_differences(input_vcf, populations, population_ids, maxmissing=None):
	positions = []
	fixed_differences = []
	with VariantFile(input_vcf) as vcf:
		for rec in vcf.fetch():
			position = rec.pos
			rec_alleles = rec.alleles
			alleles = {pop: [] for pop in population_ids}
			fixed = False
			for pop in population_ids:
				for sample in populations[pop]:
					a1 = rec.samples[sample]['GT'][0]
					a2 = rec.samples[sample]['GT'][1]
					if not a1 == None:
						alleles[pop].append(a1)
					if not a2 == None:
						alleles[pop].append(a2)
			missingness_exceeded = False
			if maxmissing is not None:
				for pop in population_ids:
					if not (len(alleles[pop]) / 2) / len(population[pop]) > maxmissing * len(population[pop]):
						missingness_exceeded = True
			else:
				for pop in population_ids:
					for pop in population_ids:
						if not len(alleles[pop]) / 2 == len(populations[pop]):
							missingness_exceeded = True
			pop1 = population_ids[0]
			pop2 = population_ids[1]
			if not missingness_exceeded:
				if len(set(alleles[pop1])) == 1 and len(set(alleles[pop2])) == 1 and not alleles[pop1][0] == alleles[pop2][0]:
					fixed_differences.append({"position": position, population_ids[0]:rec_alleles[alleles[pop1][0]],population_ids[1]:rec_alleles[alleles[pop2][0]]})
	return fixed_differences

parser = argparse.ArgumentParser(description="Extract SNPs that are different between two specified populations.")

parser.add_argument('-i', '--input-vcf', type=str, help='Input vcf (can be gzipped)', required=True)
parser.add_argument('--popfile', type=str, help="File with sample in first column and population/species assignment in second, separated by tab.")
parser.add_argument('-o', '--output', type=str, help='Output tsv file.')

args = parser.parse_args()

#parse input arguments
input_vcf = args.input_vcf
output = args.output
popfile = args.popfile

#parse popfile
population_ids, populations = parse_popfile(popfile)

#check that number of populations/species is 2:
if not len(population_ids) == 2:
	print("Number of populations/species must be two.")
	sys.exit()

#identify fixed differences
fixed_differences = find_fixed_differences(input_vcf, populations, population_ids, maxmissing=None)

#write to output
with open(output, "w") as of:
	of.write("Position\t" + population_ids[0] + "_variant\t" + population_ids[1] + "_variant\n")
	for diff in fixed_differences:
		of.write(str(diff["position"]) + "\t" + diff[population_ids[0]] + "\t" + diff[population_ids[1]] + "\n")
	of.close()
