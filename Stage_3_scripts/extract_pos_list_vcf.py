#module load pysam/0.16.0.1-python3.8.7

from pysam import VariantFile
import os
import sys

#path to vcf-file
input_vcf=str(sys.argv[1])

#output prefix, "positions.tsv" will be appended to this
output=str(sys.argv[2])

with VariantFile(input_vcf) as vcf:
    outfile=open(os.path.join(output + "_pos_list.tsv"),'w')
    outfile.write("\t".join(["pos"]) + "\n")
    for rec in vcf.fetch():
        outfile.write(str(rec.pos) + "\n")

