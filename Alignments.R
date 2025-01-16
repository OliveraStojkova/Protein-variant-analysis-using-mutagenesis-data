### Compare native DNA with the reference DNA used above

native_dna <- readDNAStringSet("native_DNA.fa")

# Pairwise alignment
# Local alignment
alignment_dna_local <- pairwiseAlignment(subject = wild_type_sequence, pattern = native_dna, type = "local")

# Global alignment
alignment_dna_global <- pairwiseAlignment(subject = wild_type_sequence, pattern = native_dna, type = "global")

# Show alignments
alignment_dna_local
alignment_dna_global

# Translate native DNA sequence
native_protein <- toString(translate((native_dna)))
native_protein_nostring <- translate(native_dna)
wild_type_protein_nostring <- translate(wild_type_sequence)

# Alignment - Protein level
# Local alignment
alignment_protein_local <- pairwiseAlignment(pattern = wild_type_protein_nostring, subject = native_protein_nostring, type = "local")

# Global alignment 
alignment_protein_global <- pairwiseAlignment(pattern = wild_type_protein_nostring, subject = native_protein_nostring, type = "global")

# Show alignments
alignment_protein_local
alignment_protein_global
