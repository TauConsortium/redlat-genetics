#!/bin/bash

## Ensure that both input filenames only include '.' at the end of the filename (ie. .vcf.gz)
## any other '.' in the filename will affect the naming in the concordance output

# check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <array_vcf>.vcf.gz <genome_vcf>.vcf.gz"
	exit 1
fi

module load cluster/R/4.1.0
module load cluster/htslib

current_dir=$(pwd)

# extract vcf basenames
array_vcf_basename=${1%.vcf.gz}
genome_vcf_basename=${2%.vcf.gz}

# unzip genome file
gunzip -k ${genome_vcf_basename}.vcf.gz

# array to store job ids
job_ids=()

# submitting filter and concordance jobs
for af in 0_01 0_001 0_0001 0_00001 0_000001 0_0000001; do
	# submit sbatch job to filter based on af
	filter_job_id=$(sbatch --parsable --job-name filter_vcf_${af} -o filter_vcf_${af}.out /cluster/home/jtaylor/scripts/Check-Concordance/filter_vcf_job.sh ${array_vcf_basename}.vcf.gz $af)

	# run snpsift concordance on each af
	concordance_job_id=$(sbatch --parsable --job-name snpsift_conc_${af} -o snpsift_concordance_${af}.out --dependency=afterok:$filter_job_id /cluster/home/jtaylor/scripts/Check-Concordance/snpsift_concordance_job.sh ${genome_vcf_basename}.vcf ${array_vcf_basename}_AF-${af}.vcf $af)

	# store concordance job id
	job_ids+=($concordance_job_id)
done

#cCreate a dependency string
dependency_string=$(IFS=,; echo "${job_ids[*]}")

# submit a dummy job that waits for all previous jobs to complete
dummy_job_id=$(sbatch --parsable --job-name waiting --dependency=afterany:$dependency_string --wrap="echo 'All jobs completed'")

# wait until dummy job completes
echo "Waiting for all jobs to complete..."
while squeue -j $dummy_job_id | grep -q $dummy_job_id; do
	sleep 300
done

echo "All SLURM jobs have completed."

# create concordance files for R script
for af in 0_01 0_001 0_0001 0_00001 0_000001 0_0000001; do
	cat concordance_${genome_vcf_basename}_${array_vcf_basename}_AF-${af}.by_sample.txt | awk '{ print $1 "\t" $14 "\t" $15 "\t" $16 "\t" $19 "\t" $20 "\t" $21 "\t" $24 "\t" $25 "\t" $26 }' > AF-${af}_conc.tsv;
done

# run R script to calculate and combine concordances
Rscript /cluster/home/jtaylor/scripts/Check-Concordance/combine_concordance.R AF-0_01_conc.tsv AF-0_001_conc.tsv AF-0_0001_conc.tsv AF-0_00001_conc.tsv AF-0_000001_conc.tsv AF-0_0000001_conc.tsv $current_dir

# remove all temp files
rm ${genome_vcf_basename}.vcf
for af in 0_01 0_001 0_0001 0_00001 0_000001 0_0000001; do
	rm ${array_vcf_basename}_AF-${af}.vcf
done
