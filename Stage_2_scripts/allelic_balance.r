AB=read.table("interval_24_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("pre_gatk_interval_2_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("post_gatk_interval_2_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_y_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("chr21_101456_ML_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("new_concatenated_gatk_filt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_gatk_filt_bcf_filt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_gatk_filt_repet_filt_bcf_filt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("AB_stats_concatenated_gatk_filt_repet_filt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_gatk_filt_repet_filt_bcf_filt_NODEPTH_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)


AB=read.table("concatenated_gatk_filt_repet_filt_bcf_filt_DEPTH2_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_gatk_filt_repet_filt_bcf_filt_DEPTH3.vcf.gz_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)


AB=read.table("concatenated_gatk_filt_repet_filt_bcf_filt_DEPTH5_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB=read.table("concatenated_nofilt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)


AB=read.table("concatenated_gatk_filt_rm_bcf1_test1_repet_filt_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)


AB=read.table("concatenated_gatk_filt_rm_repet_filt_bcf1tst2_bcf2_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) +
  facet_wrap(~sample)

AB1=read.table("concatenated_subsets_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB1, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) + xlim(0,1) +
  facet_wrap(~sample)

AB2=read.table("concatenated_subsets_preAB_allele-balance.tsv", header=T)
library(ggplot2)
ggplot(AB2, aes(x=allele_balance, group=sample)) + 
  geom_density(adjust=1.5) + xlim(0,1) +
  facet_wrap(~sample)

