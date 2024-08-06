# Check-Concordance

## Description
This pipeline is designed to check the concordance between imputed arrays and whole genome sequencing (WGS). It involves filtering the imputed array VCF file based on varying allele frequencies (AFs), comparing two VCF files (filtered imputed array and WGS) at various AF thresholds, and then using an R script to calculate and combine the concordance data.

## Requirements
- bcftools
- htslib
- Java (for SnpSift)
- R
- SLURM Workload Manager (for job memory and CPU allocation)

## Installation
Ensure that all required software (bcftools, htslib, Java, and R) is installed on your system. The script also uses specific SLURM commands for job scheduling, so SLURM needs to be available in your cluster environment.

## Usage
Run the script with two VCF files (array and genome) as arguments. Both files should be compressed (`.vcf.gz`), indexed, and should not contain periods in their filenames except for the file extension, this will affect the naming for the SnpSift output files.

```bash
./vcf_concordance_analysis.sh <array_vcf>.vcf.gz <genome_vcf>.vcf.gz
```

## Arguments
- <array_vcf>.vcf.gz - The imputed array VCF file.
- <genome_vcf>.vcf.gz - The WGS VCF file.

## Output
The script generates several intermediate files and ultimately outputs a TSV file that contains the concordance data for various allele frequencies. The final output is output_concordance.tsv which is filtered based on specific criteria defined in the R script.

## Script Components
- vcf_concordance_analysis.sh: The main script that runs the workflow.
- filter_vcf_job.sh: A SLURM job script for filtering the VCF based on allele frequency.
- snpsift_concordance_job.sh: A SLURM job script for running SnpSift concordance.
- combine_concordance.R: An R script to calculate and combine concordances.

## Additional Notes
- Ensure the paths to SLURM job scripts and R scripts are correct in your environment. 
- Update any SLURM job attributes to fit the needs of your analysis.

## References
- Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li. "Twelve years of SAMtools and BCFtools." GigaScience, Volume 10, Issue 2, February 2021, giab008, [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
- Cingolani, P., Patel, V.M., Coon, M., Nguyen, T., Land, S.J., Ruden, D.M., Lu, X. "Using Drosophila melanogaster as a model for genotoxic chemical mutational studies with a new program, SnpSift." Frontiers in Genetics, Volume 3, 2012, Frontiers Media SA

## Contact
For support or contributions, please contact me at:
 
`jtaylor[at]hudsonalpha.org`
