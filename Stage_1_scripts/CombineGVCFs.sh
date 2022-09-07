#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:10:00
#SBATCH -J CombineGVCFs
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# This script consolidates the contents of the intermediate GVCF. Script outputs one cohort gvcf file for each interval.
# The previous step, HaplotypeCaller, ran on all chromosomes and all samples, we will now combine through intervals and
# call the GenotypeGVCFs tool to obtain a final VCF file for each chromosome. NOTE! this includes non-variant sites.
# Be sure to have this in mind when handling the VCF file.

# ADDITIONAL IMPORTANT NOTE: Version GATK/4.1.4.1 will output stats for non-variant sites, but version GATK/4.2.0.0 will not.

#module load bioinfo-tools GATK/4
module load bioinfo-tools GATK/4.1.4.1

# give path to region file as cmd input (same as used before)
REGION=$(echo $(basename $1) | cut -d . -f 1)

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs
combine_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/COMBINEGVCFs
# mkdir -p "$out_dir"
gvcf_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/HAPLOTYPECALLER

# construct arguments file with one variant file per line (grep -v means exclude, $1 in awk means first column of the file)
find $gvcf_dir | grep $(basename $REGION) | grep -v ".tbi" > "$REGION"_variants.txt   # I replaced .idx with .tbi
awk ' { print "-V " $1 } ' "$REGION"_variants.txt > argfile."$(basename $REGION)".txt
rm -r -f "$REGION"_variants.txt
arg_file=argfile."$(basename $REGION)".txt

gatk --java-options "-Xmx24G" CombineGVCFs \
	-R "$ref_genome" \
	-O "$combine_dir"/cohort_interval_"$REGION".g.vcf.gz \
	-L $1 \
	--arguments_file $arg_file

gatk --java-options "-Xmx24G" GenotypeGVCFs \
	-R "$ref_genome" \
	-V "$combine_dir"/cohort_interval_"$REGION".g.vcf.gz \
	-L $1 \
	-O "$out_dir"/"$REGION".vcf.gz \
	--include-non-variant-sites

