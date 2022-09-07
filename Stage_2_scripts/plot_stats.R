#clear up environment:
rm(list = ls())

#load packages:
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(optparse)

#get command inputs
option_list <- list(make_option(c("--dir"), type = 'character'))

opt <- parse_args(OptionParser(option_list = option_list))

setwd(opt$dir)

#read in data
het <- read_tsv(list.files(getwd(), patter="\\.het$"))
site_summed_depth <- read_tsv(list.files(getwd(), patter="\\.ldepth$"))
individual_depth <- read_tsv(list.files(getwd(), patter="\\.idepth$"))
individual_missing <- read_tsv(list.files(getwd(), patter="\\.imiss$"))
genotype_depth <- read_tsv(list.files(getwd(), patter="\\.ldepth.mean$"))
genotype_missing <- read_tsv(list.files(getwd(), patter="\\.lmiss$"))
genotype_quality <- read_tsv(list.files(getwd(), patter="\\.lqual$"))


#plot heterozygosity per individual:
#add heterozygosity ratio to table:
het <- het %>% add_column(heterozygosity = (het$N_SITES - het$`O(HOM)`)/het$N_SITES)
#add individual depth for sorting:
het <- het %>% full_join(individual_depth, by = "INDV")
#sort by mean depth, largest to smalles:
het <- het %>% arrange(desc(MEAN_DEPTH))
#create ordered list for plotting:
order=het$INDV
#make plot:
hetplot <- ggplot(het, aes(x=factor(INDV, level = order), y=heterozygosity)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5)) +
  xlab("Samples sorted by read depth (high-low)") + ylab("Heterozygosity") + ggtitle("Heterozygosity")

#plot depth distribution:
het <- het %>% arrange(desc(MEAN_DEPTH))
depth_plot <- ggplot(het, aes(x=factor(INDV, level = order), y=MEAN_DEPTH)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5)) +
  xlab("") + ylab("Depth") + ggtitle("Mean genotype depth per individual")

depth_plot_average <- mean(genotype_depth$MEAN_DEPTH)
depth_plot_average <- format(round(depth_plot_average, 2), nsmall = 2)
depth_plot_sd <- sd(genotype_depth$MEAN_DEPTH)
depth_plot_sd <- format(round(depth_plot_sd, 2), nsmall = 2)
paste("Average ind. dp/gt", " (", depth_plot_average, "±", depth_plot_sd,")", sep = "")
#plot individual depth per site:
site_depth_plot <- ggplot(genotype_depth, aes(genotype_depth$MEAN_DEPTH)) + geom_density() +
  xlim(0, 75) + xlab("Mean depth") + ylab("Density") + ggtitle(paste("Average ind. dp/gt", " (", depth_plot_average, "±", depth_plot_sd,")", sep = ""))


#plot summed depth per site:
summed_depth_plot_average <- mean(site_summed_depth$SUM_DEPTH)
summed_depth_plot_average <- format(round(summed_depth_plot_average, 1), nsmall = 1)
summed_depth_plot_sd <- sd(site_summed_depth$SUM_DEPTH)
summed_depth_plot_sd <- format(round(summed_depth_plot_sd, 1), nsmall = 1)


summed_depth_plot <- ggplot(site_summed_depth, aes(site_summed_depth$SUM_DEPTH)) + geom_density() +
  xlim(400, 2000) + xlab("summed depth") + ylab("Density") + ggtitle(paste("Summed dp/gt", " (", summed_depth_plot_average, "±", summed_depth_plot_sd,")", sep = ""))


#plot site quality:
site_quality_plot <- ggplot(genotype_quality, aes(QUAL)) + geom_density() +
  xlim(0, 5000) + xlab("Quality (cutoff at 5000)") + ylab("Density") + ggtitle("Genotype quality distribution")

#plot missing genotypes per individual:
missing_ind_plot <- ggplot(individual_missing, aes(x = factor(INDV, level = order), y = F_MISS)) +
  geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust= 0.5)) +
  xlab("Samples sorted by read depth (high-low)") + ylab("Missing fraction") + ggtitle("Missing genotypes per sample")


#plot missing samples per genotype:
missing_geno_plot <- ggplot(genotype_missing, aes(F_MISS)) + geom_histogram(binwidth = 0.01, aes(y=..count../sum(..count..))) +
  xlab("Missing fraction") + ylab("Fraction of genotypes") + ggtitle("Missing fraction of samples per genotype")

output <- ggarrange(site_depth_plot, site_quality_plot, missing_geno_plot, summed_depth_plot, hetplot, depth_plot, missing_ind_plot, ncol = 4, nrow = 2, heights = c(1,2))

outfile <- paste(str_replace(list.files(getwd(), pattern ="\\.het$"), pattern =".het$", ""), ".pdf", sep = "")

setwd("..")

pdf(file = outfile, width = 18, height = 9)
output
dev.off()
