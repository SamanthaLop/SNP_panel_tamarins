#!/bin/bash -l

#SBATCH -A PROJ
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 24:00:00
#SBATCH -J vcffix
#SBATCH -e /PROJ_DIR/PATH/%x-%j.error
#SBATCH -o /PROJ_DIR/PATH/%x-%j.out
#SBATCH --mail-user=MAIL@ADDRESS
#SBATCH --mail-type=FAIL

# subsetting variants from a VCF doesn't require updating INFO fields, but removing samples does. That is what vcffix does.

module load bioinfo-tools bcftools vcflib GATK/4

# Remember to update the out_dir!
#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/simpe/x_chromosome
#out_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls
#out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/individual_species_specific_SNPs/simp/interval_vcfs/4_vcffix
out_dir=/proj/proj_name/nobackup/SAM/flanking_regions/adding_to_ind_set_feb_2022/rerun_combine_genotype_gvcfs/vcffixup_5

filename=$(echo $(basename $1)| cut -d . -f 1)

# index VCF, this is needed for vcffix to run
gatk IndexFeatureFile -I $1

# use vcffixup to update INFO fields # careful, I don't think vcffixup outputs compressed vcf files.
vcffixup $1 > "$out_dir"/"$filename"_vcffixup.vcf

gzip $1
