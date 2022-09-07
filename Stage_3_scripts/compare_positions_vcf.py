#module load pysam/0.16.0.1-python3.8.7

from pysam import VariantFile
import os
import sys

simp_pos_list=str(sys.argv[1])
lwed_pos_list=str(sys.argv[2])
output=str(sys.argv[3])

# AXEL NOTES
simp_pos = []  # will create an empty list
with open(simp_pos_list) as f: # will open the simp_pos_list file
    for line in f.readlines(): # for each line in the simp_pos_list file
        simp_pos.append(line.rstrip()) # append the content to the simp_pos LIST

lwed_pos=[]
with open(lwed_pos_list) as l:
    for line in l.readlines():
        lwed_pos.append(line.rstrip())

outfile=open(os.path.join(output + "_shared_pos_list.tsv"),'w')
for item in simp_pos:
    if item in lwed_pos:
        outfile.write(str(item) + "\n")
