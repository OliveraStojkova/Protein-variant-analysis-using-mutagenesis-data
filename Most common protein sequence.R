## Determine the most common protein sequence that is not wild type and report the mutation(s) found in this sequence
# Filter out protein sequences
non_wt_protein <- mutation_data[mutation_data$protein != wt_protein, ]

# Find the most comon protein sequence
most_common_protein <- names(sort(table(non_wt_protein$protein), decreasing = TRUE) [1])

# Compare wt and mutated protein sequence and report mutations
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

mutations <- find_mutations(wt_protein, most_common_protein)

cat("The most common non-wt protein sequence is: ", most_common_protein , "\n")
cat("Mutation(s) in the most common non-wildtype protein sequence:", paste(mutations, collapse = ", "), "\n")
