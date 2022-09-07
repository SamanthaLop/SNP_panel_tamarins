#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1-00:00:00
#SBATCH -J VariantFiltration
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# GATK recommended filtering, input the chromsome file (interval_X.bed)
# Took about 20 hours at the very most
module load bioinfo-tools picard GATK/4

ref_genome=/proj/proj_name/nobackup/SAM/references/callithrix/GCA_009663435.2_Callithrix_jacchus_cj1700_1.1_genomic.fna
genotype_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs
out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/GATK_filtered


REGION=$(echo $(basename $1) | cut -d . -f 1)

gatk --java-options "-Xmx4g" VariantFiltration \
        -R $ref_genome \
        -V "$genotype_dir"/"$REGION".vcf.gz \
        -O "$out_dir"/gatk_annotated_"$REGION".vcf.gz \
        --filter-name "QD2" --filter-expression "QD<2.0" \
        --filter-name "QUAL30" --filter-expression "QUAL<30.0" \
        --filter-name "SOR3" --filter-expression "SOR>3.0" \
        --filter-name "FS60" --filter-expression "FS>60.0" \
        --filter-name "MQ40" --filter-expression "MQ<40.0"  \
        --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum<-12.5" \
        --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum<-8.0"


# set filtered gt to no call
gatk SelectVariants \
	-V "$out_dir"/gatk_annotated_"$REGION".vcf.gz \
	-O "$out_dir"/gatk_filtered_"$REGION".vcf.gz  \
	--set-filtered-gt-to-nocall

# all filters here are using the recommended thresholds by gatk
# QD is QualByDepth. This is the QUAL score normalized by allele depth (AD) for a variant.
# QUAL is phred-scaled probability that a REF/ALT polymorphism exists at this site given sequencing data. A 10 indicated a 1 in 10 chance of error
#   and a 100 indicates a 1 in 10^10 chance. Not often very useful property for evaluating the quality of a variant call.
# SOR is the StrandsOddsRatio, a way to estimate strand bias.
# FS is FisherStrand, the Phred-scaled probability that there is strand bias at the site. (Less bias-> closer to 0)
# MQ is RMSMappingQuality, the root mean square mapping qualiy over all the reads at the site.
# MQRankSum is the MappingQualityRankSumTest, u-based z-approximation from the Rank Sum Test for mapping qualities. It compares the mapping
#    qualities of the reads supporting the reference allele and the alt allele.
# ReadPosRankSum is a u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference 
#    and alternate alleles are diff within the reads.

