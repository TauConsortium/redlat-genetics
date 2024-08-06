#!/bin/bash

# check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <vcf_array>.vcf.gz <vcf_genome>.vcf.gz"
    exit 1
fi

module load cluster/bcftools
module load cluster/htslib

job_mem=$SLURM_MEM_PER_NODE
job_mem=$((job_mem * 80 / 100))

vcf_array_basename=${1%.vcf.gz}
vcf_genome_basename=${2%.vcf.gz}

bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf_array_basename}.vcf.gz -Oz -o ${vcf_array_basename}.new_id.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ${vcf_genome_basename}.vcf.gz -Oz -o ${vcf_genome_basename}.new_id.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

tabix -p vcf ${vcf_array_basename}.new_id.vcf.gz
tabix -p vcf ${vcf_genome_basename}.new_id.vcf.gz

# intersect the two vcfs, and make a list of the variants that are in both
bcftools isec -n=2 ${vcf_array_basename}.new_id.vcf.gz ${vcf_genome_basename}.new_id.vcf.gz -o vcf.intersect --threads ${SLURM_JOB_CPUS_PER_NODE}

# convert bcftools output to a list of ids in the chrom:pos:ref:alt format
cat vcf.intersect | awk '{ print $1 ":" $2 ":" $3 ":" $4 }' > intersect.filter

# filter both vcfs for intersecting variants
bcftools view -i 'ID=@intersect.filter' ${vcf_array_basename}.new_id.vcf.gz -Oz -o ${vcf_array_basename}_intersect.new_id.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}
bcftools view -i 'ID=@intersect.filter' ${vcf_genome_basename}.new_id.vcf.gz -Oz -o ${vcf_genome_basename}_intersect.new_id.vcf.gz --threads ${SLURM_JOB_CPUS_PER_NODE}

# index new vcfs
tabix -p vcf ${vcf_array_basename}_intersect.new_id.vcf.gz
tabix -p vcf ${vcf_genome_basename}_intersect.new_id.vcf.gz

plink2 --vcf ${vcf_array_basename}_intersect.new_id.vcf.gz --make-bed --out ${vcf_array_basename} --set-missing-var-ids @:#:\$r:\$a --const-fid --new-id-max-allele-len 50 truncate --vcf-half-call missing --max-alleles 2 --memory ${job_mem}
plink2 --vcf ${vcf_genome_basename}_intersect.new_id.vcf.gz --make-bed --out ${vcf_genome_basename} --set-missing-var-ids @:#:\$r:\$a --const-fid --new-id-max-allele-len 50 truncate --vcf-half-call missing --max-alleles 2 --memory ${job_mem}

plink --bfile ${vcf_array_basename} --bmerge ${vcf_genome_basename} --make-bed --memory ${job_mem} --out redlat_merged_array_genome

plink2 --bfile redlat_merged_array_genome --keep-allele-order --recode vcf-iid --out redlat_merged_array_genome --memory ${job_mem}

bgzip -@${SLURM_JOB_CPUS_PER_NODE} redlat_merged_array_genome.vcf
tabix -p vcf redlat_merged_array_genome.vcf.gz
