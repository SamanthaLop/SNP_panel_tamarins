#KING_plots
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_csv("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/subsets/KING/simperator_king/simperator_king.kin0", sep='\t')
df = pd.read_csv("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/subsets/KING/lweddelli_king/lweddelli_king.kin0",sep="\t")


LWED = [101455.0,101456.0,101457.0,101458.0,101459.0,101460.0,101461.0,101462.0,101463.0,101465.0,101466.0,101469.0,101470.0,101471.0,101472.0,101473.0,101475.0]
SIMP = [101467.0,101468.0,101477.0,101478.0,101479.0,101480.0,101481.0]
lwed_list = []
simp_list = []

for index, row in df.iterrows():
    if (row["ID1"] in LWED) and (row["ID2"] in LWED):
        lwed_list.append(row)
    elif (row["ID1"] in SIMP) and (row["ID2"] in SIMP):
        simp_list.append(row)
        
lwed_df = pd.DataFrame(lwed_list)
simp_df = pd.DataFrame(simp_list)

df_lwed=lwed_df.pivot_table(columns="ID1", index = "ID2", values ="Kinship").reset_index()

# replace neg values with 0
df_lwed[df_lwed < 0] = 0
df_lwed_2=df_lwed.set_index("ID2")

# find max
s=df_lwed_2.select_dtypes(include=[np.number]).max()
max=s.max()

df_simp=simp_df.pivot_table(columns="ID1", index = "ID2", values ="Kinship").reset_index()

# replace neg values with 0
df_simp[df_simp < 0] = 0
df_simp_2=df_simp.set_index("ID2")

# find max
s=df_simp_2.select_dtypes(include=[np.number]).max()
max=s.max()

fig, ax = plt.subplots(1,1,figsize=(20,15))
#cmap=sns.cm.rocket_r
cmap= sns.color_palette("Blues", as_cmap=True)
heatmap=sns.heatmap(df_lwed_2, vmin=0, vmax=max, ax=ax, annot=True, cmap=cmap)
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=15)
heatmap.set_yticklabels(heatmap.get_yticklabels(),rotation =45, horizontalalignment='right', fontsize=15)
for text in ax.texts:
    text.set_size(14)
    if float(text.get_text()) >= 0.177:
        text.set_size(20)
        text.set_weight('bold')
        text.set_color("black")
ax.set_xlabel("Sample ID", fontsize=20)
ax.set_ylabel("Sample ID", fontsize= 20)
ax.collections[0].colorbar.set_label("Estimated kinship coefficient",fontsize=15)
#ax.set_title("Pairwise kinship coefficients for saddleback tamarin samples", fontsize=25)
plt.show()

fig, ax = plt.subplots(1,1,figsize=(20,15))
#cmap=sns.cm.rocket_r
cmap= sns.color_palette("Blues", as_cmap=True)
heatmap=sns.heatmap(df_simp_2, vmin=0, vmax=max, ax=ax, annot=True, cmap=cmap)
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=15)
heatmap.set_yticklabels(heatmap.get_yticklabels(),rotation =45, horizontalalignment='right', fontsize=15)
for text in ax.texts:
    text.set_size(14)
    if float(text.get_text()) >= 0.177:
        text.set_size(20)
        text.set_weight('bold')
        text.set_color("black")

ax.set_xlabel("Sample ID", fontsize=20)
ax.set_ylabel("Sample ID", fontsize= 20)
ax.collections[0].colorbar.set_label("Estimated kinship coefficient",fontsize=15)
#ax.collections[0].colorbar.set_
#ax.set_title("Pairwise kinship coefficients for emperor tamarin samples", fontsize=25)
plt.show()

simp_df
x,y = [],[]
x.append (simp_df.IBS0)
y.append (simp_df.Kinship)
fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(111)
ax.plot(x,y,'o',color='black')
ax.hlines(0.0442,0, 0.10, color='black')
ax.hlines(0.0884,0, 0.10, color='black')
ax.hlines(0.177,0, 0.10, color='black')
ax.hlines(0.354,0, 0.10, color='black')
plt.xlabel("Proportion of SNPs with zero Identical By State sharing", fontsize=20)
plt.ylabel("Estimated pairwise kinship coefficient", fontsize=20)
plt.text(0.007,0.358,"Duplicate/MZ twin", fontsize= 15)
plt.text(0.006,0.180,"1st degree relationship", fontsize= 15)
plt.text(0.006,0.091,"2nd degree relationship", fontsize= 15)
plt.text(0.006,0.047,"3rd degree relationship", fontsize= 15)
plt.show()

lwed_df
x,y = [],[]
x.append (lwed_df.IBS0)
y.append (lwed_df.Kinship)
fig = plt.figure(figsize=(20,15))
ax = fig.add_subplot(111)
ax.plot(x,y,'o',color='black')
ax.hlines(0.0442,0, 0.10, color='black')
ax.hlines(0.0884,0, 0.10, color='black')
ax.hlines(0.177,0, 0.10, color='black')
ax.hlines(0.354,0, 0.10, color='black')
plt.xlabel("Proportion of SNPs with zero Identical By State sharing", fontsize=20)
plt.ylabel("Estimated pairwise kinship coefficient", fontsize=20)
plt.text(0.006,0.358,"Duplicate/MZ twin", fontsize= 15)
plt.text(0.006,0.180,"1st degree relationship", fontsize= 15)
plt.text(0.006,0.091,"2nd degree relationship", fontsize= 15)
plt.text(0.006,0.047,"3rd degree relationship", fontsize= 15)
plt.show()
