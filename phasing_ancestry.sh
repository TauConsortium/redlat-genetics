#!/bin/bash
# Script for Phasing with SHAPEIT5 and determining Local ancestry with RFMix2
# Juliana Acosta-Uribe 2024

## run script as:
# chmod u+x phasing_ancestry.sh
# nohup ./phasing_ancestry.sh > phasing_ancestry.log


## REQUIRED SOFTWARE:
# bcftools: https://samtools.github.io/bcftools/bcftools.html#index
# shapeit5: https://odelaneau.github.io/shapeit5/
# RFMix2: https://github.com/slowkoni/rfmix

## REQUIRED FILES:
# 1.  File you want to phase (must be a file.bcf)
# You can convert vcf --> bcf using bcftools:
# `bcftools convert --output-type b file.vcf > file.bcf`
query_file='example.chr17' # Basename (prefix) for input .bcf files

# 2. Fasta file of reference genome. Must have an index (.fai) in same directory 
# fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
# The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html
fasta_file='hg38.fa.gz'

# 3. Genetic map
# Three columns: bp, Chromosome, cM
g_map='chr17.b38.gmap.gz'
# if you do not provide this file, 

# 4. Reference panel for Ancestry determination
# reference_ancestry file should be phased data. Can be a .vcf.gz or .bcf 
# Must have an index (.tbi or .csi) in same directory 
reference_ancestry='1KGP.afr-eur-eas-nat.samples_chr17.bcf'

# 5. Sample Map
# The sample map file specifies which subpopulation each reference sample represents. 
# It is tab delimited text with at least two columns. 
# The first column gives the sample name or identifier, which must match the one used in the reference VCF/BCF. 
# The second column is a string naming a subpopulation and may contain spaces (e.g., "European", or "East_African").
ancestry_sample_map='afr-eur-eas-nat.sample-map.txt'


## OPTIONAL FILES;
# 1. Reference Panel for phasing
# File should be phased data. Can be a .vcf.gz or .bcf. Must have an index (.tbi or .csi) in same directory 
reference_phasing='CCDG_14151_B01_GRM_WGS_2020-08-05_chr17.filtered.shapeit2-duohmm-phased.vcf.gz' 
# If you don't want to use a reference for phasing, delete `--reference ${reference_phasing}` from the phase_common command 

# 3. Pedigree file
# This file contains one line per sample having parent(s) in the dataset and three columns (kidID fatherID and motherID), separated by TABs for spaces.Use NAs for unknown parents (in the case of duos). In output file, the first offspring haplotype is transmitted by the father, the second by the mother.
pedigree_file='example.ped'
# If the family relationships are specified in a plink style family.fam. you can generate a pedigree_file using awk
# awk '{ if ($3 == 0) $3 = "NA"; if ($4 == 0) $4 = "NA"; print $2, $3, $4}' file.fam > pedigree_file
# If you don't want to use a pedigree for phasing, delete `--pedigree ${pedigree_file}` from the phase_common command 
## VARIABLES
chromosome='chr17' # Cromosome that is being processed
threads='24' 

## 1. PREPARE YOUR FILES

## Align to reference and normalize
# Chromosomes must be encoded as chr#. 
# If you need to change chromosome designation from # to chr# you need to create a file (e.g. chromosomes.txt) where every line has two columns: 'old name' 'new name' (e.g. 1 chr1)
#bcftools annotate --rename-chrs chr14.txt --output-type b ${query_file}.bcf > ${query_file}.chr_id.bcf
#query_file=${query_file}'.chr_id'

bcftools norm --check-ref x --fasta-ref ${fasta_file} --threads ${threads} --output-type b ${query_file}.bcf > ${query_file}.ref.bcf
# --check-ref warn (w), exclude (x), or set/fix (s)
query_file=${query_file}'.ref'  

## Annotate with AC and AN if the file doesnt have it
bcftools +fill-tags ${query_file}.bcf --output-type b --output ${query_file}.AC.bcf -- -t AC,AN 
query_file=${query_file}'.AC'

## Index your file
bcftools index ${query_file}.bcf


## 2. PHASE YOUR FILES WITH SHAPEIT5

# 1. Phase common variants (recommended thrashold of MAF >= 0.1%)
phase_common \
--input ${query_file}.bcf \
--reference ${reference_phasing} \
--region ${chromosome} \
--map ${g_map} \
--output ${query_file}.phased.bcf \
--thread ${threads} \
--pedigree ${pedigree_file} \
--filter-maf 0.01 \
--log ${query_file}.shapeit5.log

# You can use the phased data as a scaffold and incorporate the rare variants using SHAPEIT5. 
# In most cases this is not necessary for local ancestry determination.


## USE RFMIX2 TO DETERMINE LOCAL ANCESTRY

# 1. Recode the genetic map according to RFMix2 recommendations:
# The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM.
#zcat ${g_map} | awk  '{print $2, $1, $3}' > ${g_map}.rfmix2

# 1.5 Extract a Reference panel for Ancestry determination from your Reference panel for Phasing
#bcftools view --samples-file afr-eur-eas-nat.samples.txt  --threads ${threads} --output-type b ${reference_phasing} > 1KGP.afr-eur-eas-nat.samples_chr14.bcf
#bcftools index 1KGP.afr-eur-eas-nat.samples_chr14.bcf


# 2. Run RFMix2

# rfmix 
# -f BCF/VCF file with samples to analyze.
# -r BCF/VCF file with reference individuals
# -m ${ancestry_sample_map} #Reference panel sample population classification map
# -g ${g_map}.rfmix2 #Genetic map file 
# -o ${query_file}.phased.rfmix2 #Basename (prefix) for output files 
# --chromosome=${chromosome} #Execute only on specified chromosome   
# -n 5 #Terminal node size for random forest trees recommended value
# -e 5 #Maximum number of EM iterations recommended value
# --reanalyze-reference #In EM, analyze local ancestry of the reference panel and reclassify it. Recommended flag
# --random-seed=clock 
# --n-threads=${threads}

# Make sure your ${query_file}.phased.bcf is indexed
bcftools index ${query_file}.phased.bcf

# Run Rfmix2
rfmix \
-f ${query_file}.phased.bcf \
-r ${reference_ancestry} \
-m ${ancestry_sample_map} \
-g ${g_map}.rfmix2 \
-o ${query_file}.phased.rfmix2 \
--chromosome=${chromosome} \
-n 5 \
-e 5 \
--reanalyze-reference \
--random-seed=clock \
--n-threads=${threads}

# RFMix2 will generate the following output files:
# *output.msp.tsv*: The most likely assignment of subpopulations per CRF point. Produced by computing the viterbi algorithm on the CRF. The *.msp.tsv* file is condensed such that CRF windows are combined if all query samples are in the sample subpopulations for successive windows. 
# *output.fb.tsv*: The marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point. Produced by computing the forward-backward algorithm on the CRF
# *output.sis.tsv*
# *output.Q*: RFMix2 also reports global ancestry estimated based on local ancestries identified by their algorithm in their standard output file (∗results.rfmix.Q) corresponding to the .Q output files from ADMIXTURE that can be compared to ADMIXTURE generated global ancestry directly. 
