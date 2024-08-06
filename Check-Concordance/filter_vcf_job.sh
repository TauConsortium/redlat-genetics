#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

module load cluster/bcftools

af_value=$2
af=${2//_/.}

vcf_file=$1
base_name=${vcf_file%.vcf.gz}

# run bcftools view
bcftools view -i "AF>${af}" $vcf_file -Ov -o ${base_name}_AF-${af_value}.vcf --threads 8
