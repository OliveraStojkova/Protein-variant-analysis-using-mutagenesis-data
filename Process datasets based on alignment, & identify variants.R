# Load GFP datasets
brightGFP_data <- read.table("bright_GFP_beads.counts", sep = "\t", stringsAsFactors = FALSE, , col.names = c("Counts", "Sequence"))
dimGFP_data <- read.table("dim_GFP_beads.counts", sep = "\t", stringsAsFactors = FALSE, col.names = c("Counts", "Sequence"))

# Clean up
# Length filtering
brightGFP_data <- brightGFP_data[nchar(brightGFP_data$Sequence) == nchar(native_dna), ]
dimGFP_data <- dimGFP_data[nchar(dimGFP_data$Sequence) == nchar(native_dna), ]
# Remove gaphs (-)
dimGFP_data <- dimGFP_data[!grepl("-", dimGFP_data$Sequence), ]
brightGFP_data <- brightGFP_data[!grepl("-", brightGFP_data$Sequence), ]
# Remove sequences with non-DNA letters (i.e., anything other than A, T, C, G)
dimGFP_data <- dimGFP_data[!grepl("[^ATCG]", dimGFP_data$Sequence), ]
brightGFP_data <- brightGFP_data[!grepl("[^ATCG]", brightGFP_data$Sequence), ]

# Translate sequences
dimGFP_data$protein <- as.character(lapply(dimGFP_data$Sequence, function(s) toString(translate(DNAString(s)))))
brightGFP_data$protein <- as.character(lapply(brightGFP_data$Sequence, function(s) toString(translate(DNAString(s)))))

# Remove sequences with premature stop codons (those that contain "*" anywhere in the sequence except the end)
dimGFP_data <- dimGFP_data[!grepl("\\*", dimGFP_data$protein) | regexpr("\\*", dimGFP_data$protein) == nchar(dimGFP_data$protein), ]
brightGFP_data <- brightGFP_data[!grepl("\\*", brightGFP_data$protein) | regexpr("\\*", brightGFP_data$protein) == nchar(brightGFP_data$protein), ]

# Based on the alignment filter the length of the sequences (Local protein alignment)
wt <- pattern(alignment_protein_local)
native <- subject(alignment_protein_local)

# Get start and end positions
start_wt <- start(wt)
end_wt <- end(wt)

start_native <- start(native)
end_native <- end(native)

# Trim wt and native sequences to match the alignment
wt_trimmed <- toString(subseq(wild_type_protein_nostring, start_wt, end_wt))
native_trimmed <- toString(subseq(native_protein_nostring, start_native, end_native))

# Filter the datasets so that the protein sequences match in terms of position based on the alignment 
trim_protein <- function(protein_seq, start_pos, end_pos) {
  return(substr(protein_seq, start_pos, end_pos))
}

brightGFP_data$trimmed_protein <- sapply(brightGFP_data$protein, trim_protein, start_pos = start_native, end_pos = end_native)

dimGFP_data$trimmed_protein <- sapply(dimGFP_data$protein, trim_protein, start_pos = start_native, end_pos = end_native)

# Do the same for the sarkisyan dataset - use the one after filtering unique Barcode > 1, because the lowest count in bright and dim is 2
mutation_data_filtered$trimmed_protein <- sapply(mutation_data_filtered$protein, trim_protein, start_pos = start_wt, end_pos = end_wt)

# Function to compare WT and mutant proteins
find_mutations <- function(wt_seq, mut_seq) {
  wt_aa <- unlist(strsplit(wt_seq, NULL))
  mut_aa <- unlist(strsplit(mut_seq, NULL))
  
  mutations <- c()
  for (i in seq_along(wt_aa)) {
    if (wt_aa[i] != mut_aa[i]) {
      mutation <- paste0(wt_aa[i], i, mut_aa[i]) 
      mutations <- c(mutations, mutation)
    }
  }
  return(mutations)
}

# Initialize vectors to store combined mutations
bright_combined_mutations <- c()
dim_combined_mutations <- c()

# Process the dim dataset
for (i in 1:nrow(dimGFP_data)) {
  wt_protein <- native_trimmed
  mutant_protein <- dimGFP_data$trimmed_protein[i]
  mutations <- find_mutations(wt_protein, mutant_protein)
  
  # Combine mutations into a single string, separated by commas
  if (length(mutations) > 0) {
    dim_combined_mutations[i] <- paste(mutations, collapse = ", ")
  } else {
    dim_combined_mutations[i] <- NA  # Handle cases with no mutations
  }
}

# Process the bright dataset
for (i in 1:nrow(brightGFP_data)) {
  wt_protein <- native_trimmed
  mutant_protein <- brightGFP_data$trimmed_protein[i]
  mutations <- find_mutations(wt_protein, mutant_protein)
  
  # Combine mutations into a single string, separated by commas
  if (length(mutations) > 0) {
    bright_combined_mutations[i] <- paste(mutations, collapse = ", ")
  } else {
    bright_combined_mutations[i] <- NA  # Handle cases with no mutations
  }
}

# Create new data frames with the combined mutation strings
bright_mutation_results <- data.frame(trimmed_protein = brightGFP_data$trimmed_protein, Mutations = bright_combined_mutations, stringsAsFactors = FALSE) 

dim_mutation_results <- data.frame(trimmed_protein = dimGFP_data$trimmed_protein, Mutations = dim_combined_mutations, stringsAsFactors = FALSE) 

# Merge with original datset
brightGFP_data <- cbind(brightGFP_data, bright_mutation_results)
brightGFP_data <- brightGFP_data[ , !duplicated(as.list(brightGFP_data))]%>%
  separate_rows(Mutations, sep = ",")%>%
  drop_na(Mutations)

dimGFP_data <- cbind(dimGFP_data, dim_mutation_results)
dimGFP_data <- dimGFP_data[ , !duplicated(as.list(dimGFP_data))]%>%
  separate_rows(Mutations, sep = ",")%>%
  drop_na(Mutations)

# Mutations in the Sarkisyan dataset
sarkisyan_mutations <- c()

for (i in 1:nrow(mutation_data_filtered)) {
  wt_protein <- wt_trimmed
  mutant_protein <- mutation_data_filtered$trimmed_protein[i]
  mutations <- find_mutations(wt_protein, mutant_protein)
  
  # Combine mutations into a single string, separated by commas
  if (length(mutations) > 0) {
    sarkisyan_mutations[i] <- paste(mutations, collapse = ", ")
  } else {
    sarkisyan_mutations[i] <- NA  # Handle cases with no mutations
  }
}

sarkisiyan_mutation_results <- data.frame(trimmed_protein = mutation_data_filtered$trimmed_protein, Mutations = sarkisyan_mutations, stringsAsFactors = FALSE)

mutation_data_filtered <- cbind(mutation_data_filtered, sarkisiyan_mutation_results)
mutation_data_filtered <- mutation_data_filtered[ , !duplicated(as.list(mutation_data_filtered))]

# Find unique mutations for bd_merged_GFP and Sarkisyan dataset
merge_bright_dim_datasets <- merge(brightGFP_data, dimGFP_data, by = "Mutations", suffix = c("_bright", "_dim"))

unique_mutations_bd <- unique(merge_bright_dim_datasets$Mutations)

sarkisyan <- mutation_data_filtered%>%
  separate_rows(Mutations, sep = ",")%>%
  drop_na(Mutations)
unique_mutations_sar <- unique(sarkisyan$Mutations)

# Find variants present in both datasets
vars_in_common <- intersect(unique_mutations_sar, unique_mutations_bd)

# Print the results
cat("Number of Mutations Observed in Sarkasiyan GFP Dataset and Bright-Dim Datasets:", length(vars_in_common), "\n")
cat("Number of Unique Mutations in Sarkasiyan GFP Dataset:", length(unique_mutations_sar), "\n")
cat("Number of Unique Mutations in Bright-Dim GFP Dataset:", length(unique_mutations_bd), "\n")

