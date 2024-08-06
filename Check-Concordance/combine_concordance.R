# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# list of af
af_values <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001)

i <- 1
raw_concordance <- list()

# read in all raw concordance files
for (af in af_values) {
  raw_concordance[[as.character(af)]] <- read.delim(file = args[i])
  i <- i + 1
}

new_concordance <- list()

# for all afs
for (af in af_values) {
	# count all matching calls
	raw_concordance[[as.character(af)]]$correct_calls <- raw_concordance[[as.character(af)]]$REF.REF + raw_concordance[[as.character(af)]]$ALT_1.ALT_1 + raw_concordance[[as.character(af)]]$ALT_2.ALT_2
	
	# count total number of calls
	raw_concordance[[as.character(af)]]$all_calls <- rowSums(raw_concordance[[as.character(af)]][, 2:10], na.rm = TRUE)
	
	new_col_name <- paste0("af_", as.character(af), "_concordance") # column name for concordance
	
	# calculate concordance
	raw_concordance[[as.character(af)]][[new_col_name]] <- raw_concordance[[as.character(af)]]$correct_calls / raw_concordance[[as.character(af)]]$all_calls
	
	# creat new df with sample and concordance only
	new_concordance[[as.character(af)]] <- raw_concordance[[as.character(af)]][, c("sample", new_col_name)]
}

# merge all data frames in the new_dfs list by 'sample'
merged_concordance <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), new_concordance)

# Filter out rows where af_0.01_concordance is less than 0.70
filtered_concordance <- merged_concordance[merged_concordance$af_0.0001_concordance >= 0.95, ]

# output filename
output_file <- paste0(args[7], "/output_concordance.tsv")

# write output file
write.table(filtered_concordance, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
