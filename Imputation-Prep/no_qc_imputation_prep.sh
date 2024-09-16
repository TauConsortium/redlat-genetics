#!/bin/bash

input_prefix=$1
reference_flag=$2

wraynor_script=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/HRC-1000G-check-bim.pl
snp_sift=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/SnpSift.jar

### set two reference files to hg19 or hg38 depending on user's input	
if [[ $reference_flag -eq 38 ]]
then
	freeze8_1=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg38.tab.gz
	freeze8_2=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg38.tab.gz	

elif [[ $reference_flag -eq 19 ]]
then
	freeze8_1=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg19.tab.gz
	freeze8_2=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg19.tab.gz

fi

### load modules for SspSift and plink
module load cluster/plink/1.90
module load cluster/R/4.1.1
module load cluster/htslib

### Generate .frq file.
plink \
	--freq \
	--bfile ${input_prefix} \
	--out ${input_prefix} \
	--memory ${SLURM_MEM_PER_NODE}

### run the Wraynor script to compare to ALL TOPMed freeze 8 variants to be comprehensive at this stage
perl $wraynor_script \
	-b ${input_prefix}.bim \
	-f ${input_prefix}.frq \
	-r $freeze8_1 \
	-h \
	--verbose

### The excluded genotypes are flipped to rescue more variants.
plink \
	--bfile ${input_prefix} \
	--flip Exclude-${input_prefix}-HRC.txt \
	--make-bed \
	--out ${input_prefix}_FlipExclusions \
	--memory ${SLURM_MEM_PER_NODE}

input_prefix=${input_prefix}_FlipExclusions

### generate the .frq file for this set
plink \
	--freq \
	--bfile ${input_prefix} \
	--out ${input_prefix} \
	--memory ${SLURM_MEM_PER_NODE}

### compare to PASS filter variants for TOPMed to ensure highest quality to go into imputation
perl $wraynor_script \
	-b ${input_prefix}.bim \
	-f ${input_prefix}.frq \
	-r $freeze8_2 \
	-h \
	--verbose

### ensure the run_plink script is executable
chmod 755 Run-plink.sh

### run the Run-plink.sh script
sh Run-plink.sh

### rename 23 to X for input to the TOPMed server
for i in ${input_prefix}; do sed -e "s/##contig=<ID=23/##contig=<ID=X/" ${i}-updated-chr23.vcf | awk -F $'\t' 'BEGIN {OFS = FS} { gsub("23","X",$1); print $0 }' > ${i}-updated-chrX.vcf && rm ${i}-updated-chr23.vcf; done

### zip up all VCF files
for i in *.vcf; do bgzip "$i"; done

for i in *.vcf.gz; do tabix -p vcf "$i"; done

### Make a new folder to contain only the input files for imputation
mkdir ZippedVCFsForImputationInput

mv *-updated-chr*.vcf.gz* ZippedVCFsForImputationInput/