# SNP_panel_tamarins
Scripts used in the development and validation of a GT-seq SNP panel for saddleback and emperor tamarins - Degree Project in Biology to obtain the MSc degree in Evolutionary Biology at Uppsala University. Supervised by Katerina Guschanski and Mrinalini Watsa, with the help of Axel Jensen and in collaboration with Rachel Voyt and Field Projects International.

The entire pipeline is being modulated in Snakemake to increase reproducibility.

## Stage 1 Scripts
1. fastqc.sh -> Runs fastqc analysis on raw reads.
2. multiqc.sh -> Runs multiqc on output of fastqc (or qualimap) to make a report summary of all samples.
3. ReorderFasta.sh -> Reorders the raw sequencing files so that the paired end reads match to one another. 
4. RG_data_script1.sh -> Indexes and creates the dictionary of the reference genome. It converts the fastq files into sam, adds read group data, marks adapters and reverts sam back to fastq format for mapping.
5. bwa_script2.sh -> Maps reads to reference genome with bwa mem and produces coordinate sorted bam file.
6. merge_bams_script3.sh -> Merges mapped and unmapped bam files to create a single mapped bam file with read group data.
7. dedup_script4.sh -> Takes merged bam, marks duplicates and removes them. Calls on qualimap and valsam to run on each deduplicated bam file.
8. valsam.sh -> Validate sam tool reports validity of SAM or BAM files. Quality control.
9. qualimap.sh -> Creates overview of bam file in pdf format.
10. HaplotypeCaller.sh -> Uses GATKs HaplotypeCaller to call SNPs and indels, creates intermediate gvcf file. Uses "interval files": a simple bed file (tab delimited) for each chromosome containing the chromosome name, start and stop positions for that chromosome.
11. CombineGVCFs.sh -> Uses CombineGVCF and GenotypeGVCFs to go from intermediate GVCF files to VCF files per chromosome. Includes non variant sites.

## Stage 2 scripts
1. VariantFiltration.sh -> GATK recommended hard filtering of VCF files.
2. remove_sample_vcf.sh -> removes sample that was underperforming.
3. RepeatMasker_filtration_conc.sh -> Uses bcftools to exclude variants in concatenated VCF file in repetitive regions (obtained from: https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/)
4. RepeatMasker_filtration_chrom.sh -> Uses bcftools to exclude variants in chromosome VCF file in repetitive regions.
5. Subset_VCFs.sh -> Counts number of sites, calculates r needed to obtain 100,000 or 10,000 subset of vcf.
6. Concatenate_VCFs.sh -> Uses BCFtools to concatenate a list of vcf files. Output would then be used to calculate and explore statistics for further variant filtering.
7. vcftools_stats.sh -> Takes a VCF and extracts the following stats from it: Allele frequency, mean depth of coverage per individual, mean depth of coverage for each site, site quality score for each site, proportion of missing data per sample, proportion of missing data per site, heterozygosity and inbreeding coefficient per individual, depth of coverage for each genotype (large file!), and depth per site summed across all individuals.
8. plot_stats.R -> WRITTEN BY AXEL JENSEN. Generates plots for the stats in a pdf file.
9. bcf_filtration.sh -> Filters VCF with the following filters: INFO/DP<1309 && INFO/DP>213, QUAL>30, F_MISSING<0.1, -f 'PASS'.
10. second_bcf_filtration.sh -> Uses bcftools to include only bi-allelic SNPs with MAF>0.1.
11. get_AB_stats.py -> Extracts Allelic Balance data from VCF file and outputs it in the form of columns to be able to plot.
12. allelic_balance.r -> Plots Allelic Balance data - we expect a normal distribution.
13. AB_Function.py -> WRITTEN BY AXEL JENSEN. Function. Uses Pysam, takes three inputs, the input VCF, the output VCF, and the allelic balance threshold. If AB is below threshold will set genotype to no-call.
14. AB_Filter.py -> Calls the AB_Function.
15. AB_Python_send.sh -> Calls and runs AB filter script.
16. SelectVariants.sh -> Keeps only variants with a max number of no-calls of 0.
17. heterozygosity.py -> Exploring and plotting heterozygosity per genome sample and species.

## Stage 3 scripts
1. fixeddifferences.py -> WRITTEN BY AXEL JENSEN. Finds variable sites that are fixed for each species. Requires the input VCF, the output TSV, and the samples.txt file with is a tsv list of each sample and which population (species) it belongs to.
2. Python_send_fixedDifferences.sh -> Sends fixeddifferences python script to be run on each chromosome.
3. Subset_flanking.sh -> Extracts a VCF with the list of positions obtained from fixeddifferences.
4. vcftools_thin.sh -> Filters VCF with the --thin flag so no 2 sites are within X distance from each other.
5. remove_lweddelli_vcf.sh -> Removes from a VCF all samples pertaining to Leontocebus weddelli.
6. remove_simperator_vcf.sh -> Removes from VCF all samples pertaining to Saguinus imperator.
7. vcffix.sh -> Uses vcffix - from VCFlib - to correct AC, AF, NS fields of the VCF file. This needs to be done after removing samples (not variants).
9. KING.sh -> Generates PLINK binary format set of files from a VCF input file and runs KING -kinship on them.
10. remove_rel_simp.sh -> Removes SIMP related individuals from VCF file.
11. remove_rel_lwed.sh -> Removes LWED related individuals from VCF file.
12. MAF_filtering_vcf.sh -> Filters out any variants from VCF with a MAF below 0.3 or above 0.45.
13. extract_pos_list_vcf.sh -> Produces list of positions from a VCF.
14. Python_send_extract_pos_list_vcf.sh -> Sends the previous script to be run.
15. compare_positions_vcf.py -> WRITTEN BY AXEL JENSEN. Locates and produces the overlapping variants between two lists (SIMP and LWED).
16. Python_send_compare_positions_vcf.sh -> Sends the previous script to be run.
17. fixeddifferences_sex1_with_genotypes.py -> Finds sites on the X chromosome where females are homozgous and males are heterozygous. Also outputs expected genotypes for each.
18. remove_males.sh -> Removes either lwed of simp males from input VCF.
19. obtain_flanking.sh -> Uses bedtools flank to extract a list of flanking regions from a bed file list of positions.
20. consecutive_positions_in_vcf.py -> Looks for consecutive regions above a certain input threshold within the input VCF file. The script will look in the consecutive region for the SNP from an input SNP list and output the SNP, the start, stop, and length of the consecutive region.
21. add_brackets_to_fasta.py -> Adds a [] around the SNP in fasta files.
22. VCFtoMSA2_plusconsensus.py -> Modified from a script written by AXEL JENSEN. Creates a Multiple Sequence Alignment from a VCF and a consensus sequence fasta file using IUPAC codes. 
23. KING_plots.py -> Python scripts to plot king output.

## Stage 4 scripts
1. ANGSD_consensus_seq.sh -> Uses ANGSD to obtain consensus fasta from input bam file. As it stands now, it uses IUPAC codes for variable sites.
2. mfeprimer.sh -> Runs MFEprimer - a software for in silico PCR - using a fasta genome (UNZIPPED!) as template and input primers.
3. splitting_mfeoutput.sh -> Splits the output of MFEprimer into 6 files: primer list, hairpin list, dimer list, description of potential amplicons, amplicon details, and parameters/citation/website/contact/time used.
4. filter_amplicons.py -> Filters amplicon description file with default 65 Tm upper limit, 55 Tm lower limit and 1000bp as upper size limit. The script can receive inputs (-u -l -s, respectively) if these need to change though.
5. filter_dimers.py -> Filters the dimer list tsv MFE primer output.
6. initial_in_silico_pcr.py -> Script -for my reference, not exactly intended for others- interpreting amplicon results of MFE primer.
7. GT_pipeline.sh -> Sends all the necessary GTseq scripts to be run. GTseq scripts are in the GTseq_scripts directory and also available at: https://github.com/GTseq

## Stage 5 scripts
1. initial_exploration_gt_output.py -> Python scripts used for initial exploration of GT seq output - what variants/samples produced zero reads, how well the rest of them performed etc. Purely exploratory script.
2. detect_sex_from_gt_output.py -> Python script to find genotypes (male or female) for selected sexSNPs in the genotyped samples.
3. detect_spp_from_gt_output.py -> Python script to find genotypes (lwed or simp) for selected spSNPs in the genotypes samples.
4. ind_id_explore_1.py -> Python script for initial exploration of indSNPs - part 1
5. ind_id_explore_2.py -> Python script for initial exploration of indSNPs - part 2
6. twin_sex_genotypes.py -> Python script that explores possible patterns of sexSNP genotypes with respect to the sex of the twin of the focal sample.
7. hethom_homhom_distances.py -> Python script to explore distances between hom hom sites and het het sites in indSNPs.
