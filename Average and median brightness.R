### Calculate average brightness for each mutation
# determine differences to the native protein sequence. For simplicity consider each position independently, i.e.,
# regardless of whether this mutation was only observed in context of other mutations. So, if a protein with A23P and S84Q is
# reported to have a medianBrightness of 3.5, then the brightness of A23P is 3.5 and the brightness of S84Q is 3.5 (in that protein).
# Aaverage the brightness across all contexts, so if i.e. there were 3 unique(!) proteins containing the A23P mutation,
# the averaged brightness of A23P is mean(brightness(protein1),brightness(protein2),brightness(protein3))

# As a control for the averaged data across different sequences, create a subset of the dataset 
# where only single-mutation sequences are considered.

# Calculate differences between wt and mutant sequences
substitution_distances <- function(s1, s2) {
  mapply(function(c1, c2) sum(c1 != c2), strsplit(s1, ""), strsplit(s2, ""))
}

# Add number of differences from WT to data
mutation_data$substitutions <- mapply(substitution_distances, mutation_data$protein, MoreArgs = list(wt_protein))

# Function to determine the differences between wild-type and mutated protein sequences
get_mutations <- function(wt_seq, mut_seq) {
  wt_split <- strsplit(wt_seq, "")[[1]]
  mut_split <- strsplit(mut_seq, "")[[1]]
  
  mutations <- mapply(function(wt_aa, mut_aa, pos) {
    if (wt_aa != mut_aa) {
      return(paste0(wt_aa, pos, mut_aa))
    }
    return(NA)
  }, wt_aa = wt_split, mut_aa = mut_split, pos = 1:length(wt_split))
  
  return(na.omit(mutations))
}

mutations_list <- list()

# Loop over each mutated protein sequence, find mutations, and store details
for (i in seq_len(nrow(mutation_data))) {
  mut_seq <- mutation_data$protein[i]
  mut_brightness <- mutation_data$medianBrightness[i]
  
  # Find mutations compared to wild-type
  mut_info <- get_mutations(wt_protein, mut_seq)
  
  if (length(mut_info) > 0) {
    # Create a data frame for mutations and their brightness
    mut_df <- data.frame(
      mutation = mut_info,
      brightness = mut_brightness,
      stringsAsFactors = FALSE
    )
    # Add the mutation data to the list
    mutations_list[[i]] <- mut_df
  }
}

# Combine all the mutation data frames into one
mutations_df <- bind_rows(mutations_list)

# Calculate the average brightness for each mutation across different protein contexts
mutation_summary <- mutations_df %>%
  group_by(mutation) %>%
  summarise(
    avg_brightness = mean(brightness),
    std_err = sd(brightness) / sqrt(n()),
    n = n()
  )

# Subset single mutants
single_mutants <- subset(mutation_data, substitutions == 1)

# Add a mutation column to indicate which mutation it is
single_mutantions <- single_mutants %>%
  rowwise() %>%
  mutate(mutation = {
    # Compare with the wild-type sequence
    wt_seq <- wt_protein
    mut_seq <- protein

    # Find the mutation by comparing sequences
    mutation_info <- sapply(seq_along(strsplit(wt_seq, '')[[1]]), function(pos) {
      if (substring(wt_seq, pos, pos) != substring(mut_seq, pos, pos)) {
        # Get the amino acid change in A124B format
        paste0(substring(wt_seq, pos, pos), pos, substring(mut_seq, pos, pos))
      } else {
        NA
      }
    })

    # Remove NAs and return as a string
    mutation_string <- na.omit(mutation_info)
    paste(mutation_string, collapse = "; ")  
  })

single_mutants_df <- as.data.frame(single_mutantions)

### Compare the median brightness of single-mutation sequences to the averaged data

# Merge the two dataframes 
comparison_df <- single_mutants_df %>%
  merge(mutation_summary, by = "mutation") 

# Plot 
ggplot(comparison_df, aes(x = medianBrightness, y = avg_brightness)) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_brightness - stdErr, ymax = avg_brightness + stdErr), width = 0.1) +
  geom_errorbarh(aes(xmin = medianBrightness - stdErr, xmax = medianBrightness + stdErr), height = 0.1) +
  labs(
    title = "Median Brightness vs. Average Brightness",
    x = "Single Mutation Median Brightness",
    y = "Average Brightness"
  ) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  theme_minimal() 

### Compare the median brightness of single-mutation sequences to the averaged data - filtered for data with more than 1 unique barcode

# Compare the medianBrightness of those single-mutant sequences to the averaged data you created above. In the example
# above that would mean comparing the brightness of the A23P single mutant protein to the average brightness over all proteins that
# contain an A23P mutation.

# Filter the data to only keep the data for sequences that have more than 1 unique barcode
mutation_data_filtered <- mutation_data%>%
  filter(uniqueBarcodes > 1)

mutations_filtered_list <- list()

# Loop over each mutated protein sequence in the filtered dataset
for (i in seq_len(nrow(mutation_data_filtered))) {
  mut_seq <- mutation_data_filtered$protein[i]
  mut_brightness <- mutation_data_filtered$medianBrightness[i]
  
  # Find mutations compared to wild-type sequence
  mut_info <- get_mutations(wt_protein, mut_seq)
  
  # Format mutations in A124B format and store details
  if (length(mut_info) > 0) {
    # Create a data frame for mutations and their brightness
    mut_df <- data.frame(
      mutation = mut_info,
      brightness = mut_brightness,
      stringsAsFactors = FALSE
    )
    
    # Add the mutation data to the list
    mutations_filtered_list[[i]] <- mut_df
  }
}

# Combine all the mutation data frames into one
mutations_filtered_df <- bind_rows(mutations_filtered_list)

# Calculate the average brightness for each mutation across different contexts again
mutation_summary_filtered <- mutations_filtered_df %>%
  group_by(mutation) %>%
  summarise(
    avg_brightness = mean(brightness),
    std_err = sd(brightness) / sqrt(n()),
    n = n()
  )

# Subset single mutants
single_mutants_filtered <- subset(mutation_data_filtered, substitutions == 1)

# Add a mutation column to indicate which mutation it is
single_mutants_filtered <- single_mutants_filtered %>%
  rowwise() %>%
  mutate(mutation = {
    # Compare with the wild-type sequence
    wt_seq <- wt_protein
    mut_seq <- protein

    # Find the mutation by comparing sequences
    mutation_info <- sapply(seq_along(strsplit(wt_seq, '')[[1]]), function(pos) {
      if (substring(wt_seq, pos, pos) != substring(mut_seq, pos, pos)) {
        # Get the amino acid change in A124B format
        paste0(substring(wt_seq, pos, pos), pos, substring(mut_seq, pos, pos))
      } else {
        NA
      }
    })

    # Remove NAs and return as a string
    mutation_string <- na.omit(mutation_info)
    paste(mutation_string, collapse = "; ")  
  })

single_mutants_filtered_df <- as.data.frame(single_mutants_filtered)

comparison_filtered_df <- single_mutants_filtered_df %>%
  merge(mutation_summary_filtered, by = "mutation") 

# Plot
ggplot(comparison_filtered_df, aes(x = medianBrightness, y = avg_brightness)) +
  geom_point() +
  geom_errorbar(aes(ymin = avg_brightness - stdErr, ymax = avg_brightness + stdErr), width = 0.1) +
  geom_errorbarh(aes(xmin = medianBrightness - stdErr, xmax = medianBrightness + stdErr), height = 0.1)+
  labs(
    title = "Median Brightness vs. Average Brightness (Unique Barcode > 1)",
    x = "Single Mutation Median Brightness",
    y = "Average Brightness"
  ) +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  theme_minimal()
