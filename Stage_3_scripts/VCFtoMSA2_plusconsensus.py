import os,sys
from pysam import VariantFile
import argparse
import random

#input_vcf="/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/rerun_combine_genotype/one_vcf_per_snp/new_homologous_vcf_post_gatk/CM018939.1:900750-900947_snp_900829.vcf.gz"
#outfile="/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/rerun_combine_genotype/one_vcf_per_snp/trying_to_make_VCFtoMSA_output_include_hets_as_IUPAC/CM018939.1:900750-900947_snp_900829.fasta"
#consensus_fasta="/proj/proj_name/nobackup/SAM/flanking_regions/female_hom_male_het_SNPs/rerun_combine_genotype/one_vcf_per_snp/trying_to_make_VCFtoMSA_output_include_hets_as_IUPAC/CM018939.1:900750-900947_snp_900829_consensus.fasta"

# this script takes a VCF and creates an MSA and consensus sequence fasta files using the extended IUPAC codes. I allows for up to three alleles in a single position.
# HOW TO RUN: python3 <script.py> -i <input.vcf> -o <output.fasta> -con <output_consensus.fasta>


def openFile(filepath, mode='r'):
    if filepath.endswith(".gz"):
        return bgzf.open(filepath, 'wt')
    else:
        return open(filepath, 'w')

def closeFile(f, filepath):
    if filepath.endswith(".gz"):
        return bgzf.close(f)
    else:
        return f.close()

def IUPAC(alleles):
        IUPAC = {
                'M': ['A','C'],
                'R': ['A','G'],
                'W': ['A','T'],
                'S': ['C','G'],
                'Y': ['C','T'],
                'K': ['G','T'],
                'B': ['C','G','T'],
                'D': ['A','G','T'],
                'H': ['A','C','T'],
                'V': ['A','C','G']
        }
        hit = False
        for c in IUPAC:
                if set(IUPAC[c]) == set(alleles):
                        code = c
                        hit = True
        if hit == False:
                code = "N"
        return code

def get_bases(rec, samples=None, use_ambiguities=True, set_filtered_to="N", ignore_filters=False):
        import random
        sample_alleles = {}
        if samples==None:
                samples = list(rec.samples)
        alleles = {i: allele for i,allele in enumerate(rec.alleles)}
        longest_allele = 1
        for a in alleles:
                if len(alleles[a]) > longest_allele:
                        longest_allele = len(alleles[a])
        for a in alleles:
                if len(alleles[a]) < longest_allele:
                        alleles[a] = alleles[a] + "-" * (longest_allele - len(alleles[a]))
                if alleles[a] == "*" or alleles[a] == None:
                        alleles[a] = "-" * longest_allele
        if not ''.join(list(rec.filter)) == 'PASS':
                if not ignore_filters:
                        for sample in samples:
                                sample_alleles[sample] = set_filtered_to
                        return sample_alleles
        for sample in list(rec.samples):
                adip = rec.samples[sample]['GT']
                if None in adip:
                        a = "-" * longest_allele
                        sample_alleles[sample] = a
                else:
                        if len(set(adip)) > 1:
                                sample_alleles[sample] = [alleles[i] for i in adip]
                                if use_ambiguities:
                                        varying_length = False
                                        for i in sample_alleles[sample]:
                                                if len(i) > 1:
                                                        varying_length = True
                                        if varying_length:
                                                print("Warning: " + sample + " is heterozygote for variants of varying lengths at " + str(rec.pos) + ", will randomly choose one of the alleles to output. ")
                                                a = sample_alleles[sample][random.randint(0,1)]
                                        else:
                                                a = IUPAC(sample_alleles[sample])
                                        sample_alleles[sample] = a
                        else:
                                a = alleles[adip[random.randint(0,1)]]
                                sample_alleles[sample] = a
        return sample_alleles


def WriteSeqDictToFasta2(seqDict, outfile, append=False):
    with openFile(outfile) as f:
        for name,seq in seqDict.items():
            f.write('>' + name + "\n" + seq + "\n")
        closeFile(f, outfile)

# getting consensus
def ConsensusSequence(seqDict, iupac_threshold=0, missingness_threshold=1):
    #check that all sequences are of the same length
    seqlengths = set([len(seq) for seq in seqDict.values()])
    if not len(seqlengths) == 1:
        print("Must be an alignment with equal sequence lengths.")
        sys.exit()
    alnLen = list(seqlengths)[0]
    consensus = {'consensus_sequence': ""}
    stats = []
    #loop through all positions in the alignment and get the consensus
    for i in range(0,alnLen):
        bases = [seqDict[seq][i] for seq in seqDict.keys()]
        #check all alleles present and put them in a dict for counting
        alleles = {allele: 0 for allele in list(set(bases))}
        #count them
        for allele in alleles.keys():
            alleles[allele] = bases.count(allele)
        #count missingness
        missingness = sum(v for k,v in alleles.items() if k in ['N','n','-']) / len(seqDict)
        #remove gaps and Ns from dict and change counts to frequencies
        alleles = {k: v / len(seqDict) for k, v in alleles.items() if not k in ['N','n','-']}
        #sort to get most common allele first
        alleles = {k: v for k, v in sorted(alleles.items(), key=lambda item: item[1], reverse=True)}
        if len(alleles.keys()) == 1:
            consensus_base = list(alleles.keys())[0]
        else:
            if len(alleles) > 3:
                print("More than three alleles present at position " + str(i) + ". Exiting.")
                sys.exit()
            elif len(alleles) == 3 or len(alleles) == 2:
                #checking if the minor allele freq (pos 1 in dict) passes the threshold for iupac code (set to 0 by default)
                if alleles[list(alleles.keys())[1]] >= iupac_threshold:
                    consensus_base = IUPAC(list(alleles.keys()))
        if missingness > missingness_threshold:
                    consensus_base = "N"
        stats.append({'pos':str(i),'missingness': missingness, 'alleles':alleles})
        consensus['consensus_sequence'] = consensus['consensus_sequence'] + consensus_base
    return consensus, stats


#vcfFile="/Users/axeljensen/test.vcf.gz"

parser = argparse.ArgumentParser(description="Some functions for converting a vcf file to an alignment. Note: differing lengths on indels will be calibrated with '-', but those won't necessary come out aligned.")

parser.add_argument('-i', '--input', type=str, help='Input vcf', required=True)
parser.add_argument('-o', '--output', type=str, help='Output file in fasta.', required=True)
parser.add_argument('-con', '--consensus', type=str, help='Consensus output file', required=True)
parser.add_argument('-r', '--region', type=str, help='Optionally, a region can be given as chrom:start-end')

parser.add_argument('--samples', type=str, help="File with samples to include, one per line. Defaults to all samples in vcf.")

parser.add_argument('--ignore-filters', action="store_true",help="If given, sites will be output regardless of filters in the INFO-column. Default is changing such sites to 'N'")

parser.add_argument('--set-filtered-to', type=str, default="N", help="Optionally give a custom letter for changing filtered sites to")

parser.add_argument('--handle-heterozygotes', type=str, default="iupac", help="iupac or random, random will sample a random allele in case of heterozygosity, iupac uses iupac ambiguity codes.")

args = parser.parse_args()



#parse command inputs
input_vcf = args.input
outfile = args.output
consensus_fasta=args.consensus
if args.region:
        chrom,interval=args.region.split(":")[0],args.region.split(":")[1]
        start=int(interval.split("-")[0])
        end=int(interval.split("-")[1])
        region = [chrom,start,end]
else:
        region = [None,None,None]
if args.samples:
        samples = []
        with open(args.samples) as f:
                for line in f.readlines():
                        if not (line == ""):
                                samples.append(line.rstrip())
else:
        samples=None
if args.ignore_filters:
        ignore_filters = True
else:
        ignore_filters = False

set_filtered_to = args.set_filtered_to

if args.handle_heterozygotes == "iupac":
        use_ambiguities = True
else:
        use_ambiguities = False


with VariantFile(input_vcf) as vcf:
        for rec in vcf.fetch(region[0],region[1],region[2]):
                bases = get_bases(rec, samples, use_ambiguities=use_ambiguities, set_filtered_to=set_filtered_to, ignore_filters=ignore_filters)
                try:
                        for sample,base in bases.items():
                                seqDict[sample] = seqDict[sample] + base
                except NameError:
                        seqDict = {sample: base for sample,base in bases.items()}

consensus=ConsensusSequence(seqDict)[0]

WriteSeqDictToFasta2(seqDict, outfile, append=False)

WriteSeqDictToFasta2(consensus, consensus_fasta, append=False)
