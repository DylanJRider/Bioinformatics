#Dylan Rider
#Lab 10 Protein

# installing packages
install.packages("UniprotR")
install.packages("protti")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GenomicAlignments")

#load libraries 
library(Biostrings)
library(UniprotR)
library(protti)

# dna to protein
sequence_dna <- readDNAStringSet("sequence.fasta")
protein_seq <- translate(sequence_dna)

names(protein_seq) <- "sequence_protein"
writeXStringSet(protein_seq, "protein_sequence.fasta")

#accessions
accessions <- c("P0A799", "P08839")
writeLines(accessions, "uniprot_accessions.txt")

## accessions <- readLines("uniprot_accessions.txt")
## accessions <- trimws(accessions)
## accessions <- accessions[accessions != ""]

accession_string <- paste (accessions, collapse = ",")

go_info <- GetProteinGOInfo(accessions) 

#plot results
PlotGoInfo(go_info)

PlotGOAll(
  GOObj = go_info,
  Top = 10,
  directorypath = getwd(),
  width = 8,
  height = 5
)

#pathology and disease information
pathology_info <- GetPathology_Biotech(accessions)
disease_info <- Get.diseases(pathology_info)
#Structural information
write.csv(disease_info, "sequence_disease_info.csv", row.names = FALSE)
uniprot_data <- fetch_uniprot(accessions)
write.csv(uniprot_data, "sequence_uniprot_info.csv", row.names = FALSE)
# pdb 
pdb_ids <- c("1ZMR", "2HWG")
pdb_data <- fetch_pdb(pdb_ids)
write.csv(pdb_data, "sequence_pdb_data.csv", row.names = FALSE)
#alphafold 
alphafold_data <- fetch_alphafold_prediction(uniprot_ids = accessions, return_data_frame = TRUE)
write.csv(alphafold_data, "sequence_alphafold_data.csv", row.names = FALSE)


