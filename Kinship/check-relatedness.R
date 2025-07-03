# Script to check disclosed vs. genetic relatedness
# Juliana Acosta-Uribe July 2025


# 1. Setup and Data Import
progeny_file = "MA234.txt"    # pedigree data exported from progeny
king_file = "king.kin.txt"         # KING kinship output

progeny_data = read.delim(
  file        = progeny_file,
  header      = TRUE,
  sep         = "\t",
  check.names = FALSE,
  na.strings  = c("", "NA"),
  stringsAsFactors = FALSE
)

# 2. Clean Column Names
colnames(progeny_data) = gsub(" ", "_", colnames(progeny_data))
progeny_data[] = lapply(progeny_data, function(x) if (is.character(x)) gsub("\\s+", "", x) else x)

# 3. Construct Key Identifiers
progeny_data$IID = ifelse(
  !is.na(progeny_data$PIDN) & progeny_data$PIDN != "",
  progeny_data$PIDN,
  progeny_data$Individual_Name
)
progeny_data$Father_ID = paste(progeny_data$Pedigree_Name,
                               progeny_data$Father_ID, sep = "_")
progeny_data$Mother_ID = paste(progeny_data$Pedigree_Name,
                               progeny_data$Mother_ID, sep = "_")

# 4. Derive Parental PIDs
lookup = function(parent_col) {
  rel_zero = paste(progeny_data$Pedigree_Name, "0", sep = "_")
  idx      = match(progeny_data[[parent_col]], progeny_data$Individual_Name)
  pidn_val = progeny_data$PIDN[idx]
  ifelse(
    progeny_data[[parent_col]] == rel_zero, 0,
    ifelse(
      !is.na(pidn_val) & pidn_val != "",
      pidn_val,
      progeny_data[[parent_col]]
    )
  )
}
progeny_data$PID = lookup("Father_ID")
progeny_data$MID = lookup("Mother_ID")

# 5. Encode Sex
progeny_data$Sex = ifelse(
  progeny_data$Sex == "M", 1,
  ifelse(progeny_data$Sex == "F", 2, 3)
)

# 6. Build Pedigree and Kinship
if (!requireNamespace("kinship2", quietly = TRUE)) install.packages("kinship2", quiet = TRUE)
library(kinship2)

progeny_fixed = fixParents(
  id = progeny_data$IID,
  dadid = progeny_data$PID,
  momid = progeny_data$MID,
  sex = progeny_data$Sex,
  missid = 0
)
ped = with(progeny_fixed, pedigree(id = id, dadid = dadid, momid = momid, sex = sex, missid = "0"))

# Generate a kinship matrix based on disclosed relationships
kinship_matrix = kinship(ped)

# Plot the pedigree
plot(ped, cex=0.2)

# 7. Load KING Output and Compare
king = read.delim(
  file = king_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ids_king = sort(unique(c(king$ID1, king$ID2)))
missing_ids = setdiff(ids_king, progeny_fixed$id)
if (length(missing_ids) > 0) {
  warning("Individuals in KING data not in pedigree: ", paste(missing_ids, collapse = ", "))
} else {
  message("All KING IDs are represented in pedigree")
}
idx1 = match(king$ID1, rownames(kinship_matrix))
idx2 = match(king$ID2, colnames(kinship_matrix))
king$kinship2 = NA_real_
valid = !is.na(idx1) & !is.na(idx2)
king$kinship2[valid] = kinship_matrix[cbind(idx1[valid], idx2[valid])]

# 8. Identify parent-offspring duos in KING data
po_idx = mapply(
  function(x, y) {
    cond1 = x %in% progeny_fixed$id &&
      (y == progeny_fixed$dadid[progeny_fixed$id == x] || y == progeny_fixed$momid[progeny_fixed$id == x])
    cond2 = y %in% progeny_fixed$id &&
      (x == progeny_fixed$dadid[progeny_fixed$id == y] || x == progeny_fixed$momid[progeny_fixed$id == y])
    cond1 | cond2
  },
  king$ID1, king$ID2
)
po_duos = data.frame(ID1 = king$ID1[po_idx], ID2 = king$ID2[po_idx])
message("Expected parent-offspring duos in KING data:")
print(nrow(po_duos))
print(po_duos)

# 9. Identify full-sibling duos in KING data
fs_idx = mapply(
  function(x, y) {
    idx_x = which(progeny_fixed$id == x)
    idx_y = which(progeny_fixed$id == y)
    length(idx_x) == 1 & length(idx_y) == 1 &
      progeny_fixed$dadid[idx_x] == progeny_fixed$dadid[idx_y] &
      progeny_fixed$momid[idx_x] == progeny_fixed$momid[idx_y]
  },
  king$ID1, king$ID2
)
fs_duos = data.frame(ID1 = king$ID1[fs_idx], ID2 = king$ID2[fs_idx])
message("Expected full-sibling duos in KING data:")
print(nrow(fs_duos))
print(fs_duos)

# 10. Identify 2nd degree duos with kinship2 in [0.0884, 0.177]
# Only include pairs where both IDs are in ids_king
sd_idx = which(
  kinship_matrix >= 0.0884 & kinship_matrix <= 0.177,
  arr.ind = TRUE
)
sd_idx = sd_idx[
  rownames(kinship_matrix)[sd_idx[,1]] %in% ids_king &
    colnames(kinship_matrix)[sd_idx[,2]] %in% ids_king,
  , drop = FALSE
]
sd_duos = data.frame(
  ID1      = rownames(kinship_matrix)[sd_idx[,1]],
  ID2      = colnames(kinship_matrix)[sd_idx[,2]],
  kinship2 = kinship_matrix[sd_idx]
)
message("Expected second degree duos in KING data:")
print(nrow(sd_duos))

# 11. Assign 'pedigree_rel' based on detected duos and kinship thresholds
king$pedigree_rel = ifelse(
  paste(king$ID1, king$ID2) %in% paste(po_duos$ID1, po_duos$ID2) |
    paste(king$ID2, king$ID1) %in% paste(po_duos$ID1, po_duos$ID2),
  "PO",
  ifelse(
    paste(king$ID1, king$ID2) %in% paste(fs_duos$ID1, fs_duos$ID2) |
      paste(king$ID2, king$ID1) %in% paste(fs_duos$ID1, fs_duos$ID2),
    "FS",
    ifelse(
      paste(king$ID1, king$ID2) %in% paste(sd_duos$ID1, sd_duos$ID2) |
        paste(king$ID2, king$ID1) %in% paste(sd_duos$ID1, sd_duos$ID2),
      "2nd",
      ifelse(
        king$kinship2 >= 0.0442 & king$kinship2 < 0.0884,
        "3rd",
        "UN"
      )
    )
  )
)

# 12. Consistency Check
# Add a column where you check if the disclosed and genotypic relationships are the same
king$rel_check = king$InfType == king$pedigree_rel

# 13. Add 'warnings' column: flag mismatches where expected relationship is PO, 1st, or 2nd
king$warnings = ifelse(
  !king$rel_check & (king$InfType %in% c("PO", "1st", "2nd") |
                       king$pedigree_rel %in% c("PO", "1st", "2nd")),
  "warning",
  NA_character_
)

# 14. Preview Results
head(king)

# 15. Export updated KING data as a tab-delimited file
write.table(
  king,
  file = "king.kin.Kinship2-updated.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 16. Generate 'expected_found' summary dataframe and export
expected_found = data.frame(
  id = ids_king,
  Expected_PO  = sapply(ids_king, function(id) sum(po_duos$ID1 == id | po_duos$ID2 == id)),
  Found_PO     = sapply(ids_king, function(id) sum((king$ID1 == id | king$ID2 == id) & king$InfType == "PO")),
  Expected_FS  = sapply(ids_king, function(id) sum(fs_duos$ID1 == id | fs_duos$ID2 == id)),
  Found_FS     = sapply(ids_king, function(id) sum((king$ID1 == id | king$ID2 == id) & king$InfType == "FS")),
  Expected_2nd = sapply(ids_king, function(id) sum(sd_duos$ID1 == id | sd_duos$ID2 == id)),
  Found_2nd    = sapply(ids_king, function(id) sum((king$ID1 == id | king$ID2 == id) & king$InfType == "2nd"))
)
head(expected_found)

write.table(
  expected_found,
  file = "disclosed_vs_genetic_kinship.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 17. Plot all pair-wise kinship differences Pedigree Kinship - Genetic Kinship

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", quiet = TRUE)
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2", quiet = TRUE)
library(ggplot2)
library(reshape2)

# Generate a matrix with differences
diff_mat = matrix(NA, nrow = length(ids_king), ncol = length(ids_king),
                  dimnames = list(ids_king, ids_king))

for (i in seq_len(nrow(king))) {
  id1 = king$ID1[i]
  id2 = king$ID2[i]
  if (id1 %in% ids_king & id2 %in% ids_king) {
    diff_val = king$kinship2[i] - king$Kinship[i]
    diff_mat[id1, id2] = diff_val
    diff_mat[id2, id1] = diff_val
  }
}
# Keep only the lower triangle
diff_mat[upper.tri(diff_mat)] = NA

diff_df = na.omit(melt(diff_mat, varnames = c("ID1", "ID2"), value.name = "diff"))

# Plot using ggplot2
p = ggplot(diff_df, aes(x = ID1, y = ID2, fill = diff)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Kinnship values, Pedigree (kinship2) - Genetic (KING)",
    x = "", y = "", fill = "Difference"
  )

plot(p)

# Save the heatmap as a PNG file
ggsave(
  filename = "kinship_difference_heatmap.png",
  plot = p,
  width = 8,
  height = 6,
  units = "in"
)
