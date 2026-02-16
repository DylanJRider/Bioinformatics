#packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("msa")

install.packages("seqinr")
install.packages("phangorn")

library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)

#Reading FASTA 
seq1 <- readDNAStringSet("egfr_flank.fasta")
seq2 <- readDNAStringSet("sequence-2.fasta")
seq3 <- readDNAStringSet("sequence-3.fasta")
seq4 <- readDNAStringSet("sequence-4.fasta")
seq5 <- readDNAStringSet("sequence-5.fasta")
seq6 <- readDNAStringSet("sequence-6.fasta")
sequences <- c(seq1, seq2, seq3, seq4, seq5, seq6)
sequences
#combine sequences
names (sequences) <- c("EGFR_flank",
                       "Acomys_spinosissimus_Z96068",
                       "Eubalaena_glacialis",
                       "Melontaenia_splendida",
                       "Anguilla_rostrata",
                       "Acomys_spinosissimus_Z96067")
print(sequences)

#Muscle alignment
alignment <- msa(sequences, method = "Muscle")

print(alignment, show = "complete")

aligned <- as(alignment, "DNAStringSet")

#Gaps
gap_total <- sum(letterFrequency(aligned, "-"))
gap_total

#Alignment length
align_length <- width(as(alignment, "DNAStringSet"))[1]
align_length

#GC Content
base_counts <- colSums(alphabetFrequency(aligned)[, c("A", "T", "G", "C")])
GC_percent <- (base_counts["G"] + base_counts["C"]) /
  sum(base_counts) * 100
GC_percent

#Convert to seqinr 
alignment_seqinr <- msaConvert(alignment, type = "seqinr::alignment")

#distance matrix
dist_matrix <- as.matrix(
  dist.alignment(alignment_seqinr, "identity")
  )
dist_matrix

#Translate first sequence to amino acids
original_seq <- sequences[["Acomys_spinosissimus_Z96068"]]
trimmed_seq <- subseq(original_seq, start = 1, width = floor(width(original_seq)/3)*3)
AA_sequence <- Biostrings::translate(trimmed_seq, genetic.code = Biostrings::getGeneticCode("SGC1"))
AA_sequence

#Export alignment to FASTA
alignment_phy <- msaConvert(alignment, type= "phangorn::phyDat")

write.phyDat(alignment_phy, file = "cytb_alignment_output.fasta", format = "fasta")

