### Find variants after matching the protein with the local protein alignment 

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

# Filter sarkisyan dataset so that only sequences with mutations that are present in both datasets are shown
mutation_data_expanded <- mutation_data_filtered %>%
  separate_rows(Mutations, sep = ",")%>%
  drop_na(Mutations)

common_mutations_df <- mutation_data_expanded %>%
  filter(Mutations %in% vars_in_common)

common_mutations_df <- common_mutations_df %>%
  group_by(Mutations) %>%
  summarize(no_observed = n(),
            avg_brightness = mean(medianBrightness, na.rm=TRUE)) %>%
  ungroup()

## Plot average brightness vs. log(brigh/dim) for variants present in both datasets
# Calculate log(bright/dim) ratio

merged_data_bright_dim_filtered <- merge_bright_dim_datasets %>%
  filter(Mutations %in% vars_in_common)

merged_data_bright_dim_filtered <- merged_data_bright_dim_filtered %>%
  group_by(Mutations) %>%
  summarize(Count_in_Bright = mean(Counts_bright, na.rm = TRUE),
            Count_in_Dim = mean(Counts_dim, na.rm = TRUE)) %>%
  ungroup()

# Pseudocount value (to avoid division by zero)
pseudocount <- 1

# Total counts for bright and dim
total_bright <- sum(merged_data_bright_dim_filtered$Count_in_Bright) + pseudocount * nrow(merged_data_bright_dim_filtered)

total_dim <- sum(merged_data_bright_dim_filtered$Count_in_Dim) + pseudocount * nrow(merged_data_bright_dim_filtered)

# Add pseudocounts to both bright and dim counts to avoid zero counts
merged_data_bright_dim_filtered$Count_in_Bright <- merged_data_bright_dim_filtered$Count_in_Bright + pseudocount
merged_data_bright_dim_filtered$Count_in_Dim <- merged_data_bright_dim_filtered$Count_in_Dim + pseudocount

# Calculate odds for each sequence using pseudocounts
merged_data_bright_dim_filtered$odds_bright <- merged_data_bright_dim_filtered$Count_in_Bright / (total_bright - merged_data_bright_dim_filtered$Count_in_Bright)

merged_data_bright_dim_filtered$odds_dim <- merged_data_bright_dim_filtered$Count_in_Dim / (total_dim - merged_data_bright_dim_filtered$Count_in_Dim)

# Calculate log-odds ratio
merged_data_bright_dim_filtered$log_odds_ratio <- log(merged_data_bright_dim_filtered$odds_bright / merged_data_bright_dim_filtered$odds_dim)

# Merge datasets to plot
merged_df <- merged_data_bright_dim_filtered %>%
  left_join(common_mutations_df %>% select(Mutations, avg_brightness), by = "Mutations")

# Plot
ggplot(merged_df, aes(y = avg_brightness, x = log_odds_ratio)) +
  geom_point(alpha = 0.7, color = "navy", size = 1) +
  labs(
    title = "Averaged Median Brightness vs. Log Odds Ratio",
    y = "Average Median Brightness",
    x = "Log Odds Ratio"
  ) +
  theme_minimal()
