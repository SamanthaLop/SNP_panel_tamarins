from AB_Function import min_AB_filter
import sys

# calls upon min_AB_filter, a function contained in another python script, goal is to set genotypes to no call if below 0.25 allelic balance.
# threshold can be changed as you see fit

# (input_vcf, output_vcf, AB_threashold)
min_AB_filter(sys.argv[1], sys.argv[2], 0.25)
