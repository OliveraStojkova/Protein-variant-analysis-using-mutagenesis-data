### Incorporating MAVE assay data and GnomAD data 

abundance_data <- read.table("prism_mave_036_VKOR_abundance.txt", header = TRUE)
activity_data <- read.table("prism_mave_035_VKOR_ab_activity.txt", header = TRUE)
gnomad_data <- read.csv("gnomad_new.csv")

# Scatterplot of all variants described by the two MAVEs
abundance_activity_merged_data <- merge(abundance_data, activity_data, by = "variant", suffixes = c("_abundance", "_activity"))

ggplot(abundance_activity_merged_data, aes(x = score_abundance, y = score_activity))+
  geom_point()+
  labs(
    title = "Abundance vs. Activity of Variants Described by MAVEs",
    x = "Abundance Score",
    y = "Activity Score"
  ) +
  theme_minimal()
  
# Scatterplot of only variants listed in gnomAD                                           
# Function to convert 3-letter amino acid code to 1-letter code
three_to_one <- function(aa) {
  aa_map <- c(
    "Ala" = "A", "Arg" = "R", "Asn" = "N", "Asp" = "D", "Cys" = "C",
    "Gln" = "Q", "Glu" = "E", "Gly" = "G", "His" = "H", "Ile" = "I",
    "Leu" = "L", "Lys" = "K", "Met" = "M", "Phe" = "F", "Pro" = "P",
    "Ser" = "S", "Thr" = "T", "Trp" = "W", "Tyr" = "Y", "Val" = "V",
    "Ter" = "*"  # Stop codon
  )
  return(aa_map[aa])
}

# Process the Protein.Consequence column
gnomad_data <- gnomad_data %>%
  mutate(
    # Remove the "p." prefix
    Protein_Consequence_Cleaned = str_remove(Protein.Consequence, "^p\\."),
    
    # Extract the reference, position, and alternate amino acids
    Reference_aa = three_to_one(str_sub(Protein_Consequence_Cleaned, 1, 3)),  # Convert reference AA to 1-letter code
    Position = as.integer(str_extract(Protein_Consequence_Cleaned, "\\d+")),  # Extract position
    Alternate_aa = three_to_one(str_sub(Protein_Consequence_Cleaned, -3, -1)),  # Convert alternate AA to 1-letter code
    
    # Create the Variant column in the format <Reference><Position><Alternate>
    variant = paste0(Reference_aa, Position, Alternate_aa)  # Combine to form the variant
  )

# Filter for population variants from merged_data
population_variants <- abundance_activity_merged_data %>%
  filter(variant %in% gnomad_data$variant)

# Merge with gnomad_data to include ClinVar.Clinical.Significance
population_variants <- population_variants %>%
  left_join(gnomad_data %>% select(variant, ClinVar.Clinical.Significance), by = "variant")


# Plot score abundance vs. score phosphatase and color by clinical significance
ggplot(abundance_activity_merged_data, aes(x = score_abundance, y = score_activity)) +
  geom_point(color = "lightgray") +  # Background points for all merged_data
  geom_point(data = population_variants, aes(color = ClinVar.Clinical.Significance), size = 1) +  # Population variants colored
  geom_text_repel(data = population_variants,
                  aes(label = ifelse(ClinVar.Clinical.Significance != "Uncertain Significance", variant, '')), 
                  size = 2) +
  labs(title = "Abundance Score vs Activity Score",
       x = "Abundance Score",
       y = "Activity Score",
       color = "Clinical Significance") +
  theme_minimal() +
  theme(legend.position = "right")  


# Find intersection of the MAVEs with gnomAD
# Variants with abundance < 0.5
ab_less_than_0.5 <- population_variants %>%
  filter(score_abundance < 0.5)

# Variants with activity < 0.4
ac_less_than_0.4 <- population_variants %>%
  filter(score_activity < 0.4)

intersection_vars <- intersect(ab_less_than_0.5, ac_less_than_0.4)

# Print results
cat("Number of variants that have an abundance score less than 0.5: ", length(ab_less_than_0.5$variant), "\n")
cat("Number of variants that have an activity score less than 0.4: ", length(ac_less_than_0.4$variant), "\n")
cat("Number of variants that have an abundance score less than 0.5 and activity score less than 0.4: ", length(intersection_vars$variant), "\n")

# Variants with uncertain significance
uncretain_significance_vars <- gnomad_data%>%
  filter(ClinVar.Clinical.Significance == "Uncertain significance")

merged_uncertain_significance <- merge(uncretain_significance_vars, abundance_activity_merged_data, by = "variant")

cat("The variant labeled with uncertain significance has abundance score of: ", merged_uncertain_significance$score_abundance, "and activity score of: ", merged_uncertain_significance$score_activity, "\n")

# Benign variants
benign_variants <- gnomad_data%>%
  filter(ClinVar.Clinical.Significance == "Benign")

benign_variants$variant

# Pathogenic variants
pathogenic_variants <- gnomad_data%>%
  filter(ClinVar.Clinical.Significance == "Pathogenic")

pathogenic_variants$variant
