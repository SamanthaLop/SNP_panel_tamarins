from pysam import VariantFile
import os

#path to vcf-file
#input_vcf=

#filtered vcf output
#output=

#allele-balance threashold (specified in minor allele support)
#AB_threashold=0.25


#funtion will iterate through heterozygote positions, and if minor allele support is below threashold, the genotype will be set to nocall
def min_AB_filter(input_vcf, output_vcf, AB_threashold):
    with VariantFile(input_vcf) as vcf:
        output = VariantFile(output_vcf, 'w', header=vcf.header)
        samples = list(vcf.header.samples)
        for rec in vcf.fetch():
            if len(rec.alleles) <= 2:
                for sample in samples:
                    if rec.samples[sample]['GT'] == (0,1):
                        AB_ref,AB_alt = rec.samples[sample]['AD'][0], rec.samples[sample]['AD'][1]
                        if AB_ref <= AB_alt:
                            try:
                                AB_minor = AB_ref / (AB_ref + AB_alt)
                            except:
                                rec.samples[sample]['GT'] = (None,None)
                        else:
                            AB_minor = AB_alt / (AB_ref + AB_alt)
                        if not AB_minor >= AB_threashold:
                            rec.samples[sample]['GT'] = (None,None)
                output.write(rec)

