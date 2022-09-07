#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J bcf_filter
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script uses bcftools to further filter the variants.

module load bioinfo-tools bcftools

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/bcf1_filtered

filename=$(echo $(basename $1) | cut -d . -f 1)

#DEPTH5
#bcftools view -i 'INFO/DP<1309 && INFO/DP>213 && QUAL>30 && F_MISSING<0.1' -f 'PASS' \
# -Oz -o "$out_dir"/"$filename"_bcf.vcf.gz $1

#test1 after removing 101474 I recalculated the upper and lower DP and additionally increased the lower limit to 300 (12x per sample)
#bcftools view -i 'INFO/DP<1299 && INFO/DP>300 && QUAL>30 && F_MISSING<0.1' -f 'PASS' \
# -Oz -o "$out_dir"/"$filename"_bcf1_test1.vcf.gz $1

#test1 after removing 101474 I recalculated the upper and lower DP and additionally increased the lower limit to 300 (12x per sample)
#bcftools view -i 'INFO/DP<1299 && INFO/DP>330 && QUAL>30 && F_MISSING<0.1' -f 'PASS' \
# -Oz -o "$out_dir"/"$filename"_bcf1tst1.vcf.gz $1

#test2 increase 100 and decrease 100 from min and max limits
bcftools view -i 'INFO/DP<1299 && INFO/DP>206 && QUAL>30 && F_MISSING<0.1' -f 'PASS' \
 -Oz -o "$out_dir"/"$filename"_bcf1.vcf.gz $1

# -i means that a filter expression is coming
# -O means the output type is coming (z means compressed VCF)
# -o means that the outfile name is coming.


# Axel: simple example of filtering with bcftools. This line will keep all sites with a summed coverage across all samples < 2000 AND > 500. If the
# read coverage is too high, it can mean spurious mapping, if it is too low, it can mean unreliable variants. If the expression would
# instead have been FORMAT/DP, it would be sample-specific, but here it's INFO/DP, so it applies to the whole dataset. We also only keep SNPs
# and are restricting to sites with two alleles present (which we might need to reconsider if trying to look for SNPs that are multiallelic between
# populations but biallelic within populations.

# Axel: CAUTION. This will only output sites that pass the filters, the others will be completely exluded, this will likely be appropriate for our
# analyses but might be something to discuss down the line in case we for example want to have nonvariant sites to analyze.

# See manual page for more explanation on filter syntax: http://samtools.github.io/bcftools/bcftools.html
# $2 is the input vcf that needs filteringm

