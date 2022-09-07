#heterozygosity
import pandas as pd
import glob
from pandas import DataFrame
from statannot import add_stat_annotation
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import scipy.stats as stats



dict_of_chroms = {}

for file in glob.glob("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/interval_*_het.het"):
    df = pd.read_csv(file, sep = "\t")
    chrom = '_'.join(file.split('/')[-1].split('_')[1:2])
    dict_of_inds = {}
    for index, row in df.iterrows():
        ind = str(row["INDV"])
        obs_hom = row["O(HOM)"]
        exp_hom = row["E(HOM)"]
        n_sites = row["N_SITES"]
        f_inbreeding = row["F"]
        obs_het_sites = n_sites - obs_hom
        obs_het = obs_het_sites / n_sites
        dict_of_inds[ind] = obs_het
    dict_of_chroms[chrom] = dict_of_inds

ind_dir = DataFrame(dict_of_chroms).transpose().to_dict()

new_dict= {}
for ind in ind_dir:
    list_of_hets = []
    for chrom in ind_dir[ind]:
        list_of_hets.append(ind_dir[ind][chrom])
    new_dict[ind] = list_of_hets

SIMP_list = ["101467.0", "101468.0", "101477.0", "101478.0", "101479.0", "101480.0", "101481.0"]
LWED_list = ["101455.0", "101456.0",  "101457.0", "101458.0", "101459.0", "101460.0","101461.0", "101462.0", "101463.0", "101465.0", "101466.0","101469.0",\
     "101470.0", "101471.0", "101472.0", "101473.0",  "101475.0", "101476.0"]


SIMP_hets = []
LWED_hets = []
for ind in new_dict:
    if ind in SIMP_list:
        for value in new_dict[ind]:
            SIMP_hets.append(value)
    elif ind in LWED_list:
        for value in new_dict[ind]:
            LWED_hets.append(value)

simp_list = ["simp"]*154
lwed_list = ["lwed"]*396

species = lwed_list + simp_list
heterozygosities = LWED_hets + SIMP_hets
heterozygosities_dict = {"species": species, "het": heterozygosities}
df = pd.DataFrame(heterozygosities_dict)
simps = df.query("species == 'simp'")["het"]
lweds = df.query("species == 'lwed'")["het"]
df.groupby("het").describe()

stats.shapiro(simps)
#ShapiroResult(statistic=0.9888185858726501, pvalue=0.25829076766967773) # NORM
stats.shapiro(lweds)
#ShapiroResult(statistic=0.9962162375450134, pvalue=0.46984097361564636) # NORM
stats.levene(simps, lweds)
#LeveneResult(statistic=40.20339266282377, pvalue=4.787922173706227e-10) # DIFF

def welch_ttest(x, y): 
    ## Welch-Satterthwaite Degrees of Freedom ##
    dof = (x.var()/x.size + y.var()/y.size)**2 / ((x.var()/x.size)**2 / (x.size-1) + (y.var()/y.size)**2 / (y.size-1))
    t, p = stats.ttest_ind(x, y, equal_var = False)
    print("\n",
          f"Welch's t-test= {t:.4f}", "\n",
          f"p-value = {p:.4f}", "\n",
          f"Welch-Satterthwaite Degrees of Freedom= {dof:.4f}")

# Welch's t-test= 28.9865 
# p-value = 0.0000 
# Welch-Satterthwaite Degrees of Freedom= 206.4999




ax = sns.boxplot(x="species", y = "het", data = df)
add_stat_annotation(ax, data = df, x = "species", y = "het",  box_pairs = [("lwed", "simp")], test = "t-test_welch", text_format= "simple", loc="inside", verbose=2)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/heterozygosities_spp.png", dpi=300)



# __________
dict_of_chroms = {}

for file in glob.glob("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/interval_*_het.het"):
    df = pd.read_csv(file, sep = "\t")
    chrom = '_'.join(file.split('/')[-1].split('_')[1:2])
    dict_of_inds = {}
    for index, row in df.iterrows():
        ind = str(row["INDV"])
        obs_hom = row["O(HOM)"]
        exp_hom = row["E(HOM)"]
        n_sites = row["N_SITES"]
        f_inbreeding = row["F"]
        obs_het_sites = n_sites - obs_hom
        obs_het = obs_het_sites / n_sites
        dict_of_inds[ind] = obs_het_sites
    dict_of_chroms[chrom] = dict_of_inds

ind_dir = DataFrame(dict_of_chroms).transpose().to_dict()

new_dict= {}
for ind in ind_dir:
    list_of_hets = []
    for chrom in ind_dir[ind]:
        list_of_hets.append(ind_dir[ind][chrom])
    new_dict[ind] = list_of_hets


sum_het_sites = {}
for ind in new_dict:
    sum= 0
    for value in new_dict[ind]:
        sum = sum + value
    sum_het_sites[ind] = sum

list_of_inds = []
list_of_sum = []
list_of_sp = []
for ind in sum_het_sites:
    list_of_inds.append(ind)
    list_of_sum.append(sum_het_sites[ind])
    if ind in SIMP_list:
        list_of_sp.append("simp")
    elif ind in LWED_list:
        list_of_sp.append("lwed")


sum_dict = {"sample": list_of_inds, "species": list_of_sp, "het_sites": list_of_sum}
df = pd.DataFrame(sum_dict)

ax = sns.barplot(x = "sample", y = "het_sites", hue = "species", data = df.sort_values(by="het_sites")) #ax=ax1)
ax.tick_params(axis="x", labelrotation = 45)
ax.figure.tight_layout()
