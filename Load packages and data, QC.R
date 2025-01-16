library(ggplot2)
library(Biostrings)
library(pheatmap)
library(tidyverse)
library(reshape2)
library(tidyr)
library(dplyr)
library(seqinr)
library(pwalign)
library(ggrepel)

# Read data
mutation_data <- read.csv("nt_sequences_to_brightness.csv")
wild_type_sequence <- readDNAStringSet("avGFP_reference_sequence.fa")

######## Task 1
## Quality control
# Remove sequences that are too long, too short, or have gaps 
# Keep sequences that are the same length as the wild-type DNA sequence
mutation_data <- mutation_data[nchar(mutation_data$sequence) == nchar(wild_type_sequence), ]

# Remove gaps (-)
mutation_data <- mutation_data[!grepl("-", mutation_data$sequence), ]

## Translate sequences to protein
wt_protein <- toString(translate(wild_type_sequence))

mutation_data$protein <- as.character(lapply(mutation_data$sequence, function(s) toString(translate(DNAString(s)))))

# Remove sequences with premature stop codons (those that contain "*" anywhere in the sequence except the end)
mutation_data <- mutation_data[!grepl("\\*", mutation_data$protein) | regexpr("\\*", mutation_data$protein) == nchar(mutation_data$protein), ]

## Number of uniques barcodes, DNA and protein sequences
unique_barcodes <- length(unique(mutation_data$uniqueBarcodes))
unique_DNA_variants <- length(unique(mutation_data$sequence))
unique_protein_sequences <- length(unique(mutation_data$protein))

cat("Unique barcodes:", unique_barcodes, "\n")
cat("Unique DNA variants:", unique_DNA_variants, "\n")
cat("Unique protein sequences:", unique_protein_sequences, "\n")
