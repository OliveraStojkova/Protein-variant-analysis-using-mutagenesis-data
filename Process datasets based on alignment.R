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
