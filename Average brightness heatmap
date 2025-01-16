### Make a heatmap of average brigthness of mutant and wt amino acids

mutation_data1 <- mutation_data %>%
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

# Convert to a data frame
mutation_data1_df <- as.data.frame(mutation_data1)

# Split mutations into separate rows
# Create a new data frame with separate rows for each mutation
mutations_separated <- mutation_data1_df %>%
  separate_rows(mutation, sep = "; ") %>%  # Split by "; " to create new rows
  mutate(
    # Extract wild-type and mutant amino acids from the mutation string
    wt_aa = substring(mutation, 1, 1),  # First character is the wild-type amino acid
    position = as.numeric(substring(mutation, 2, nchar(mutation)-1)),  # Position is between 2nd and second last character
    mut_aa = substring(mutation, nchar(mutation), nchar(mutation))  # Last character is the mutant amino acid
  ) %>%
  filter(substitutions>0)


# Calculate average brightness for each position
brightness_summary <- mutations_separated%>%
  group_by(wt_aa, mut_aa) %>%
  summarise(avg_brightness = mean(medianBrightness, na.rm = TRUE), .groups = "drop")

# Create 20x20 matrix

aa_list <-c("G", "A", "V", "L", "I", "M", "F", "W", "P","S", "T", "C", "Y", "N", "Q","D", "E","K", "R", "H")

brightness_matrix <- matrix(0, nrow = length(aa_list), ncol = length(aa_list), 
                             dimnames = list(aa_list, aa_list))

# Fill in the matrix
for (i in seq_along(aa_list)) {
  for (j in seq_along(aa_list)) {
    # Extract the brightness for the specific mutation
    mutation <- paste0(aa_list[i], "->", aa_list[j])
    avg_bright <- brightness_summary %>%
      filter(wt_aa == aa_list[i], mut_aa == aa_list[j]) %>%
      pull(avg_brightness)
    
    # Fill in the matrix
    brightness_matrix[i, j] <- ifelse(length(avg_bright) > 0, avg_bright, NA)
  }
}

# Plot
brightness_df <- as.data.frame(as.table(brightness_matrix)) 
colnames(brightness_df) <- c("wt_aa", "mut_aa", "avg_brightness")

ggplot(brightness_df, aes(y = mut_aa, x = wt_aa, fill = avg_brightness)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "purple4", na.value = "lightgray") +
  labs(title = "Average Brightness Heatmap", y = "Mutant Amino Acid", x = "Wild-Type Amino Acid") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
