from pysam import VariantFile
import os

#path to vcf-file
input_vcf="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/bcf2_filtered/subsets/concatenated_subsets_preAB.vcf.gz"

#output prefix, "_allele-balance.tsv" will be appended to this
output="/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/bcf2_filtered/subsets/concatenated_subsets_preAB"

#function that will iterate through heterozygot sites and output allelebalance on a per-sample, per-site basis
def outputAlleleBalance(input_vcf, output):
	with VariantFile(input_vcf) as vcf:
		samples = list(vcf.header.samples)
		outfile = open(os.path.join(output + "_allele-balance.tsv"), "w")
		outfile.write('\t'.join(['pos','sample','allele_balance']) + "\n")
		for rec in vcf.fetch():
			for sample in samples:
				if rec.samples[sample]['GT'] == (0,1):
					AB_ref = rec.samples[sample]['AD'][0]
					AB_alt = rec.samples[sample]['AD'][1]
					try:
						AB = AB_alt / (AB_ref + AB_alt)
					except:
						AB = 0
					outfile.write(str(rec.pos) + "\t" + sample + "\t" + str(AB) + '\n')

outputAlleleBalance(input_vcf, output)
