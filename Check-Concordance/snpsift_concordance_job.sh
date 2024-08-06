#!/bin/bash

#SBATCH --partition=normal
#SBATCH --mem=120G

module load cluster/java

genome_vcf=$1
array_vcf=$2
af_value=$3


# run snpsift concordance
java -Xmx120G -jar /cluster/home/jtaylor/software/snpEff/SnpSift.jar concordance $genome_vcf $array_vcf > snpsift_conc_AF-${af_value}.txt
