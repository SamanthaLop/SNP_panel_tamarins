#!/bin/bash -l

#SBATCH -A snic2021-5-477
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J splitting_MFEoutput
#SBATCH -e /proj/proj_name/nobackup/SAM/error/%x-%j.error
#SBATCH -o /proj/proj_name/nobackup/SAM/out/%x-%j.out
#SBATCH --mail-user=samantha.lopezclinton.1915@student.uu.se
#SBATCH --mail-type=FAIL

# script that uses awk and sed to split the output of MFEprimer into 6 files for exploration

MFE_output=$1
MFE_output_filename=$(echo $(basename "$MFE_output")| cut -d . -f 1) # this might change! depends on how the mfeoutput file is called!

output_dir=/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/ANGSD_consensus_seqs/REAL_TEST/101479_round_4.2

# extract first section: PRIMER LIST:
awk '{if (match($0,"Hairpin List")) exit; print}' $MFE_output > "$output_dir"/"$MFE_output_filename"_primerlist_1.txt
awk 'NR>5 {print}' "$output_dir"/"$MFE_output_filename"_primerlist_1.txt | \
    awk 'BEGIN{print "PrimerID""\t""Seq_5_3""\t""Length_bp""\t""GC_%""\t""Tm_C""\t""Dg_kcal_mol""\t""Binding_plus""\t""Binding_minus"} ;\
    {print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' > "$output_dir"/"$MFE_output_filename"_primerlist_1.tsv

# extract second section: HAIRPIN LIST:
awk '/^Hairpin List/,/^Dimer List/' $MFE_output > "$output_dir"/"$MFE_output_filename"_pre_hairpinlist_2.txt
# remove last line
sed '$d' "$output_dir"/"$MFE_output_filename"_pre_hairpinlist_2.txt > "$output_dir"/"$MFE_output_filename"_hairpinlist_2.txt
# remove intermediate file
rm -f "$output_dir"/"$MFE_output_filename"_pre_hairpinlist_2.txt

# extract third section: DIMER LIST:
awk '/^Dimer List/,/^Descriptions of/' $MFE_output > "$output_dir"/"$MFE_output_filename"_pre_dimerlist_3.txt
# remove last line
sed '$d' "$output_dir"/"$MFE_output_filename"_pre_dimerlist_3.txt > "$output_dir"/"$MFE_output_filename"_dimerlist_3.txt
# remove intermediate file
rm -f "$output_dir"/"$MFE_output_filename"_pre_dimerlist_3.txt
awk NF "$output_dir"/"$MFE_output_filename"_dimerlist_3.txt | awk '(NR>1)' | paste -d ' ' - - - - - | sed 's/://g' | sed 's/,//g' | \
    awk 'BEGIN{print "Dimer_ID""\t""Primer_1""\t""Primer_2""\t""Score""\t""Tm_C""\t""Delta_G_kcal_mol"} ; {print$2"\t"$3"\t"$5"\t"$7"\t"$10"\t"$15}'> \
    "$output_dir"/"$MFE_output_filename"_dimerlist_3.tsv

# extract fourth section: DESCRIPTION OF POTENTIAL AMPLICONS
awk '/^Descriptions of/,/^Amplicon details/' $MFE_output > "$output_dir"/"$MFE_output_filename"_pre_descriptionamplicons_4.txt
# remove last line
sed '$d' "$output_dir"/"$MFE_output_filename"_pre_descriptionamplicons_4.txt > "$output_dir"/"$MFE_output_filename"_descriptionamplicons_4.txt
# remove intermediate file
rm -f "$output_dir"/"$MFE_output_filename"_pre_descriptionamplicons_4.txt
awk NF "$output_dir"/"$MFE_output_filename"_descriptionamplicons_4.txt | awk '(NR>3)' | awk \
    'BEGIN{print "Amp_ID""\t""HitID""\t""Size_bp""\t""FpTm_C""\t""RpTm_C""\t""FpDg_kcal_mol""\t""RpDg_kcal_mol"} ; \
    {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > "$output_dir"/"$MFE_output_filename"_descriptionamplicons_4.tsv

# extract fifth section: AMPLICON DETAILS
awk '/^Amplicon details/,/^Parameters/' $MFE_output > "$output_dir"/"$MFE_output_filename"_pre_amplicondetails_5.txt
# remove last line
sed '$d' "$output_dir"/"$MFE_output_filename"_pre_amplicondetails_5.txt > "$output_dir"/"$MFE_output_filename"_amplicondetails_5.txt
# remove intermediate file
rm -f "$output_dir"/"$MFE_output_filename"_pre_amplicondetails_5.txt
awk NF "$output_dir"/"$MFE_output_filename"_amplicondetails_5.txt | awk '(NR>1)' | awk '{if($0 ~ /^Amp/){if(s){print s;s=$0}else{s=$0}}else{s=s" "$0}}END{print s}' | \
    sed 's/-/ /g' | sed 's/:/ /g' | sed 's/:/ /g'| awk 'BEGIN{print "Amp_ID""\t""P_1""\t""P_2""\t""Chrom""\t""Start_pos""\t""Stop_pos"} ;\
    {print$2"\t"$3"\t"$5"\t"$7"\t"$8"\t"$9}' > "$output_dir"/"$MFE_output_filename"_amplicondetails_5.tsv

# merge amplicon files:
paste "$output_dir"/"$MFE_output_filename"_descriptionamplicons_4.tsv "$output_dir"/"$MFE_output_filename"_amplicondetails_5.tsv > "$output_dir"/"$MFE_output_filename"_amplicon_merged_4_5.tsv

# extract sixth section: PARAMETERS, CITATION, WEBSITE, CONTACT, TIME USED
awk '/^Parameters/,/^Placeholder/' $MFE_output > "$output_dir"/"$MFE_output_filename"_pre_parameters_6.txt
# remove last line
sed '$d' "$output_dir"/"$MFE_output_filename"_pre_parameters_6.txt > "$output_dir"/"$MFE_output_filename"_parameters_6.txt
# remove intermediate file
rm -f "$output_dir"/"$MFE_output_filename"_pre_parameters_6.txt
