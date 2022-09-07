# Twin_sex_genotypes_plotting
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from statannot import add_stat_annotation



all_samples_dict = {"1":  {"198": {"X": 41, "Y": 0}, "200": {"X": 44, "Y": 0}, "203": {"X": 60, "Y": 0}, "205": {"X": 28, "Y": 0},  "209": {"X": 126, "Y": 1},  "215": {"X": 71, "Y": 0}, "218": {"X": 68, "Y": 0}, "222": {"X": 43, "Y": 0}}, \
    "16": {"198": {"X": 5, "Y": 15}, "200": {"X": 20, "Y": 12}, "203": {"X": 16, "Y": 16}, "205": {"X": 14, "Y": 6}, "209": {"X": 52, "Y": 12}, "215": {"X": 38, "Y": 16}, "218": {"X": 4, "Y": 12}, "222": {"X":11, "Y": 1}}, \
        "8": {"195": {"X": 46, "Y": 0}, "200": {"X": 62, "Y": 0}, "205": {"X": 28, "Y": 0}, "208": {"X": 42, "Y": 0}, "209": {"X": 152, "Y": 0}, "215": {"X": 77, "Y": 0}, "222": {"X": 135, "Y":0}}, \
            "18": {"198": {"X": 70, "Y": 0}, "200": {"X": 65, "Y": 0}, "203": {"X": 65, "Y": 0}, "205": {"X": 26, "Y": 0}, "209": {"X": 135, "Y": 0}, "215": {"X": 61, "Y": 0}, "218": {"X": 50, "Y": 0}, "222": {"X": 20, "Y": 0}}, \
                "32": {"195": {"X": 71, "Y": 66}, "200": {"X": 99, "Y": 109}, "205": {"X": 89, "Y": 132}, "208": {"X": 72, "Y": 98}, "209": {"X": 199, "Y": 253}, "211": {"X": 16, "Y": 12}, "215": {"X": 62, "Y": 100}, "222": {"X": 48, "Y": 53 }}, \
                    "27": {"195": {"X": 70, "Y": 50}, "200": {"X": 117, "Y": 119}, "205": {"X": 64 , "Y": 98}, "208": {"X": 45, "Y": 61}, "209": {"X": 241, "Y": 269}, "215": {"X": 20, "Y": 65}, "222": {"X": 90, "Y": 66}}, \
                        "26": {"195": {"X": 41, "Y": 37}, "200": {"X": 90, "Y": 64}, "205": {"X": 42, "Y": 70}, "208": {"X": 44, "Y": 44}, "209": {"X": 125, "Y": 129}, "215": {"X": 28,"Y": 54}, "222": {"X": 61, "Y": 52}}, \
                            "17": {"198": {"X": 5, "Y": 12}, "200": {"X": 60, "Y": 25}, "203": {"X": 20, "Y": 26}, "205": {"X": 20, "Y": 2}, "209": {"X": 66, "Y": 20}, "215": {"X": 91, "Y": 56}, "218": {"X": 9, "Y": 20}, "222": {"X": 19,"Y": 12}}, \
                                "31": {"195": {"X": 81, "Y": 0}, "200": {"X": 103, "Y": 0}, "205": {"X": 81, "Y": 2}, "208": {"X": 78, "Y": 0}, "209": {"X": 199, "Y": 0}, "215": {"X": 68, "Y": 0}, "222": {"X": 54, "Y": 0}}, \
                                    "2": {"200": {"X": 35, "Y": 14}, "205": {"X": 18, "Y": 8}, "208": {"X": 4, "Y": 12}, "209": {"X": 44, "Y": 14}, "215": {"X": 38, "Y": 0}, "222": {"X": 14, "Y": 4}}, \
                                        "20": {"198": {"X": 15, "Y": 10}, "200": {"X": 32, "Y": 15}, "203": {"X": 14, "Y": 24}, "205": {"X": 10, "Y": 8}, "209": {"X": 56, "Y": 25}, "215": {"X": 78, "Y": 33}, "218": {"X": 14, "Y":24}, "222": {"X":15, "Y": 8}}, \
                                            "29": {"198": {"X": 248, "Y": 0}, "200": {"X": 255, "Y": 0}, "203": {"X": 159, "Y": 1}, "205": {"X": 134, "Y": 4}, "209": {"X": 395, "Y": 4}, "215": {"X": 146, "Y": 2}, "218": {"X": 162, "Y": 0}, "222": {"X": 90, "Y": 0}}, \
                                                "10": {"195": {"X": 42, "Y": 0}, "200": {"X": 66, "Y": 0}, "205": {"X": 37, "Y": 0}, "208": {"X": 40, "Y": 0}, "209": {"X": 150, "Y": 0}, "215": {"X": 67, "Y": 0}, "222": {"X": 53, "Y": 0}}, \
                                                    "23": {"198": {"X": 9, "Y": 17}, "200": {"X": 26, "Y": 12}, "203": {"X": 12, "Y": 20}, "205": {"X": 8, "Y": 6}, "209": {"X": 53, "Y": 32}, "215": {"X": 49, "Y": 10}, "218": {"X": 10, "Y": 3}}, \
                                                        "15": {"198": {"X": 24, "Y": 0}, "200": {"X": 30, "Y": 0}, "203": {"X": 32, "Y": 0}, "209": {"X": 64, "Y": 0}, "215": {"X": 37, "Y": 0}, "218": {"X": 12, "Y": 0}, "222": {"X": 17, "Y": 0}}, \
                                                            "3": {"195": {"X": 26, "Y": 0}, "200": {"X": 20, "Y": 0}, "205": {"X": 15, "Y": 0}, "208": {"X": 14, "Y": 0}, "209": {"X": 71, "Y": 0}, "215": {"X": 27, "Y": 0}, "222": {"X": 14, "Y": 0}}, \
                                                                "33": {"195": {"X": 74,"Y": 60}, "200": {"X": 86, "Y": 98}, "205": {"X": 52, "Y": 118}, "208": {"X": 62,"Y": 78}, "209": {"X": 202, "Y": 160}, "211": {"X": 4, "Y": 8}, "215": {"X": 35, "Y": 44}, "222": {"X": 40, "Y": 39}}, \
                                                                    "6": {"198": {"X": 32, "Y": 0}, "200": {"X": 54, "Y": 0}, "203": {"X": 38, "Y": 0}, "205": {"X": 20, "Y": 0}, "209": {"X": 100, "Y": 0}, "215": {"X": 75, "Y": 0}, "218": {"X": 32, "Y": 0}, "222": {"X": 32, "Y":0}}, \
                                                                        "12": {"195": {"X": 8,"Y": 10}, "200": {"X": 10, "Y": 2}, "208": {"X": 6,"Y": 16}, "209": {"X": 44, "Y": 44}, "215": {"X": 10, "Y": 20}, "222": {"X": 9, "Y": 16}}, \
                                                                            "19": {"198": {"X": 12,"Y": 4}, "200": {"X": 12, "Y": 6}, "203": {"X": 16, "Y": 20}, "205": {"X": 8, "Y": 4}, "209": {"X": 34,"Y": 24}, "218": {"X": 14, "Y": 18}, "222": {"X": 6,"Y": 7}}, \
                                                                                "34": {"197": {"X": 17, "Y": 0}, "198": {"X": 242, "Y": 0}, "200": {"X": 269, "Y": 0}, "203": {"X": 177, "Y": 0}, "205": {"X": 201, "Y": 0}, "209": {"X": 491, "Y": 0}, "215": {"X": 138, "Y": 0}, "218": {"X": 234, "Y": 0}, "222": {"X": 83, "Y": 0}}, \
                                                                                    "30": {"195": {"X": 128,"Y": 29}, "200": {"X": 212,"Y": 62}, "205": {"X": 130, "Y": 52}, "209": {"X": 325,"Y": 75}, "211": {"X": 12, "Y": 0}}, \
                                                                                        "24": {"200": {"X": 7, "Y": 6}, "203": {"X": 7,"Y": 6}, "209": {"X": 16, "Y": 14}, "218": {"X": 8, "Y": 7}, "222": {"X": 7, "Y": 3}}, \
                                                                                            "21": {"198": {"X": 45, "Y": 32}, "200": {"X": 60, "Y": 50}, "203": {"X": 41, "Y": 45}, "205": {"X":16, "Y": 22}, "209": {"X": 99, "Y": 126}, "215": {"X": 235, "Y": 56}, "218": {"X": 25, "Y": 34}, "222": {"X": 16, "Y": 13}}, \
                                                                                                "13": {"195": {"X": 12, "Y": 0}, "200": {"X": 13, "Y": 6}, "205": {"X": 4, "Y": 6}, "209": {"X": 42, "Y": 20}, "215": {"X": 4, "Y": 14}, "222": {"X": 12, "Y": 4}}, \
                                                                                                    "5": {"195": {"X": 12, "Y": 16}, "200": {"X": 24, "Y": 26}, "205": {"X": 12, "Y": 10}, "208": {"X": 18, "Y": 10}, "215": {"X": 16, "Y": 16}, "222": {"X": 14, "Y": 6}}, \
                                                                                                        "14": {"195": {"X": 24, "Y": 8}, "200": {"X": 22, "Y": 32}, "205": {"X": 18, "Y": 20}, "208": {"X": 18, "Y": 30}, "209": {"X": 58, "Y": 52}, "215": {"X": 34, "Y": 58}, "222": {"X": 27, "Y": 20}}, \
                                                                                                            "28": {"198": {"X": 193, "Y": 1}, "200": {"X": 290, "Y": 0}, "203": {"X": 278, "Y": 0}, "205": {"X": 180, "Y": 0}, "209": {"X": 560, "Y": 0}, "215": {"X": 122, "Y": 0}, "218": {"X": 149, "Y": 0}, "222": {"X": 167, "Y": 0}}, \
                                                                                                                }



females = ['8', '18', '10', '15', '3', '6', '34', '16','1','31', '29', "28", "17"]
males = ['32', '27', '26', '2', '23', '33', '19', '13', '5', '14', "12","24", "21", "30", "20"]

samples = []
chroms = []
snps = []
reads = []
sexes = []
avg_diff = []
for sample in all_samples_dict:
    for snp in all_samples_dict[sample]:
        for chrom in all_samples_dict[sample][snp]:
            read_value_X = all_samples_dict[sample][snp]["X"]
            read_value_Y = all_samples_dict[sample][snp]["Y"]
            samples.append(sample)
            snps.append(snp)
            if chrom == "X":
                chroms.append("X")
                reads.append(read_value_X)
            elif chrom == "Y":
                chroms.append("Y")
                reads.append(read_value_Y)
            if sample in males:
                sexes.append("male")
                if read_value_Y < read_value_X:
                    avg_diff.append(read_value_X - read_value_Y)
                elif read_value_Y > read_value_X:
                    avg_diff.append(read_value_Y - read_value_X)
                elif read_value_X == read_value_Y:
                    avg_diff.append(0)
            elif sample in females:
                sexes.append("female")
                if read_value_Y < read_value_X:
                    avg_diff.append(read_value_X - read_value_Y)
                elif read_value_Y > read_value_X:
                    avg_diff.append(read_value_Y - read_value_X)
                elif read_value_X == read_value_Y:
                    avg_diff.append(0)

reads_dict = {"sample": samples, "sex": sexes ,"snp": snps, "chrom": chroms, "diff_x_y": avg_diff, "reads": reads}
df = pd.DataFrame(reads_dict)
hue_order = ["X", "Y"]
sns.set(rc = {"figure.figsize":(16,6)})
ax = sns.barplot(x = "sample", y = "reads", hue = "chrom", data = df.sort_values(by="diff_x_y", ascending=False), hue_order=hue_order)
plt.ylabel("Genotyped reads")
plt.xlabel("Samples")
plt.title("Barplot of X and Y genotyped reads per sample")
ax.figure.tight_layout()
plt.legend(title="Chromosome", loc="upper right")
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/barplot_all_X_and_Y.png", dpi=300)

# now let's do only males and then only females!
females_only_df = df[df["sex"] == "female"]
x = "sample"
y = "reads"
hue = "chrom"
hue_order = ["X", "Y"]
sns.barplot(x = x, y = y, hue = hue, data = females_only_df.sort_values(by="diff_x_y", ascending=False), hue_order=hue_order)
plt.ylabel("Genotyped reads")
plt.xlabel("Samples")
plt.title("Barplot of X and Y genotyped reads per sample - females")
ax.figure.tight_layout()
plt.legend(title="Chromosome", loc="upper right")


#males
males_only_df = df[df["sex"] == "male"]
x = "sample"
y = "reads"
hue = "chrom"
hue_order = ["X", "Y"]

# Plot each sample, fems to the left and males to the right
fig, (ax1, ax2) = plt.subplots(ncols=2)
sns.barplot(x = x, y = y, hue = hue, data = females_only_df.sort_values(by="diff_x_y", ascending=False), hue_order=hue_order, ax=ax1)
ax1.set_ylabel("Average genotyped reads per sexSNP")
ax1.set_xlabel("Female samples")
ax1.set(ylim=(0,360))
ax1.legend(title="Chromosome", loc="upper right")
ax1.tick_params(axis="both", which="major", labelsize = 12)
sns.barplot(x = x, y = y, hue = hue, data = males_only_df.sort_values(by="diff_x_y", ascending=False), hue_order=hue_order, ax=ax2)
ax2.set(yticklabels=[])
ax2.set_ylabel("")
ax2.set_xlabel("Male samples")
ax2.set(ylim=[0,360])
plt.subplots_adjust(wspace=0.1)
plt.suptitle("Barplot of X and Y genotyped reads per sample")
ax2.legend(title="Chromosome", loc="upper right")
ax2.tick_params(axis="both", which="major", labelsize = 12)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/barplot_all_X_and_Y.png", dpi=300)

# okay, now let's look at FM, FF, MF, and MM plus the rest.

males = ['32', '27', '26', '23', '33', '19', '13', '5', '14']
females = ['16','1', '8', '18', '31', '29', '10', '15', '3']
females_with_fem_twin = ["28", "6", "34"]
males_with_male_twin = ["12","24", "21"]
males_with_female_twin = ["30", "20", "2"]
females_with_male_twin = ["17"]


samples = []
chroms = []
snps = []
reads = []
sexes = []

for sample in all_samples_dict:
    for snp in all_samples_dict[sample]:
        for chrom in all_samples_dict[sample][snp]:
            if chrom == "X":
                read_value = all_samples_dict[sample][snp][chrom]
                samples.append(sample)
                snps.append(snp)
                reads.append(read_value)
                chroms.append("X")
            if chrom == "Y":
                read_value = all_samples_dict[sample][snp][chrom]
                samples.append(sample)
                snps.append(snp)
                reads.append(read_value)
                chroms.append("Y")
            if sample in males:
                sexes.append("M(?)")
            elif sample in females:
                sexes.append("F(?)")
            elif sample in females_with_fem_twin:
                sexes.append("FF")
            elif sample in females_with_male_twin:
                sexes.append("FM")
            elif sample in males_with_female_twin:
                sexes.append("MF")
            elif sample in males_with_male_twin:
                sexes.append("MM")

reads_dict = {"sample": samples, "sex": sexes ,"snp": snps, "chrom": chroms, "reads": reads}

plt.clf()
df = pd.DataFrame(reads_dict)
sex_categories = ["FF","F(?)", "FM", "MF", "MM", "M(?)"]
df["sex"] = pd.Categorical(df["sex"], categories= sex_categories)
sorted_df = df.sort_values(by="sex")
hue_order = ["X", "Y"]
ax = sns.barplot(x = "sex", y = "reads", hue = "chrom", data = sorted_df, hue_order=hue_order, ci = None, estimator=sum)
ax.set_ylabel("Genotyped reads")
ax.set_xlabel("Sex categories")
ax.set_title("X and Y barplots for each sex category")
ax.legend(title="Chromosome", loc="upper right")
ax.tick_params(axis="both", which="major", labelsize = 15)

plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/barplot_sex_categories_X_and_Y.png", dpi=300)


#box plots
plt.clf()
df = pd.DataFrame(reads_dict)
sex_categories = ["FF","F(?)", "FM", "MF", "MM", "M(?)"]
df["sex"] = pd.Categorical(df["sex"], categories= sex_categories)
sorted_df = df.sort_values(by="sex")
hue_order = ["X", "Y"]
ax = sns.boxplot(x = "sex", y = "reads", hue = "chrom", data = sorted_df, hue_order=hue_order)
ax.set_ylabel("Genotyped reads per sexing SNP")
ax.set_xlabel("Sex categories")
ax.set_title("X and Y barplots for each sex category")
ax.legend(title="Chromosome", loc="upper right")
ax.tick_params(axis="both", which="major", labelsize = 15)
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/boxplot_sex_categories_X_and_Y.png", dpi=300)


