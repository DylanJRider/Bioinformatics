#package installation
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("msa")
install.packages("seqinr")
install.packages("phangorn")
 
#Load libraries

library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)

#Reading FASTA 
twenty_sequences <- readDNAStringSet("sequences.fasta")
twenty_sequences


##1 Import and align DNA sequences
#Alignment process using Muscle in the msa package
alignment <- msa(twenty_sequences, method = "Muscle")
alignment

##2 How good is the alignment? 
#The alignment can be measured by Percent Identity in order to ensure that it is a good alignment.  

#Percent Identiy 
alignment_seq <- msaConvert(alignment, type="seqinr::alignment")
d <- dist.alignment(alignment_seq, "identity")
d

#Percent Identity measures how similar two sequences are. 
#A high percent Identity index will indicate they're closely related samples, while a lower indicates low similarity. 
#A good alignment means similar sequence alignment in a way that is maximizing matches and minimizzes gaps. 


##3 Calculate Consensus sequence

#Consensus Sequence
consensusSequence <- msaConsensusSequence(alignment)
consensusSequence
consensus_dna <- DNAStringSet(consensusSequence)
writeXStringSet(consensus_dna , "consensus.fasta")

##4
# GC Content for overall alignment 
aligned <- as(alignment, "DNAStringSet")
base_counts <- colSums(alphabetFrequency(aligned)[, c("A","T", "G", "C")])
#Isolate Guanine and Cytosine 
GC_percent <- (base_counts["G"] + base_counts["C"]) / 
  sum(base_counts) * 100
GC_percent

#GC Percent is 51.56944%

##5
#Difference and find the highest difference sample:

#Convert distance object to a matrix

Dif <- as.matrix(d)

#removes self-comparisons

diag(Dif) <- NA

#Mean distance

mean_distance <- rowMeans(Dif, na.rm = TRUE)

#Highest average distance sample

highest_difference_sample <- names(which.max(mean_distance))
highest_difference_sample

#Converting alignment to a matrix format 

alignment_set <- as(alignment, "DNAStringSet")
aln_mat <-as.matrix(alignment_set)

#Most different and reference sample

ref_sample <- names(which.min(mean_distance))

sample_row <-aln_mat[highest_difference_sample, ]
ref_row <- aln_mat[ref_sample, ]

# Total differences

total_differences <- sum(sample_row != ref_row)
total_differences

#Gap(insertion/deletions) and SNPs(base substitution)

num_gaps <- sum(sample_row == "-" | ref_row == "-")
num_snps <- sum(sample_row != ref_row & 
                  !(sample_row == "-" | ref_row == "-"))

num_snps
num_gaps

#The most different individual differs from majority of samples by 
#7-8 total positions and (6-7 SNPs and 1 gap) depending on reference used

##6
#I used the concensus.fasta file and located the match.
#The match was Homo sapiens HBB gene for beta globin
#There was 0/642 0% gaps so it shows that the match is 100%.
# There was also a E-value of 0.0 and has a 100% percent Identity

##7

#Get the most different individual's Original DNA (not the aligned)
most_diff_dna <- twenty_sequences[highest_difference_sample] 

#Get sequence length 
seq_len <- width(most_diff_dna)[1]

#trim to multiple of 3
trim_length <- seq_len - (seq_len %% 3)

most_diff_dna_trimmed <- subseq(most_diff_dna, start = 1, end=trim_length)


#translate directly
protein_translation <- Biostrings::translate(most_diff_dna_trimmed)
protein_translation

#write protein to FASTA
protein_set <- AAStringSet(protein_translation)
names(protein_set) <- paste0(highest_difference_sample, "_protein")
writeXStringSet(protein_set, "most_different_protein.fasta")

##8 Protein BLAST

#The translated protein was placed into BLASTP.
#The top match was hemoglobin subunit beta HBB 
#It had a 100% identity with the lowest E-value 
#The accession number is KAI2558340.1

##9 Diseases associated with this gene
#After using the OMIM database, I discovered two main diseases:
#Sickle Cell disease
#Beta-thalassemia
#DNA codon 7 must be checked for
#sickle. Sickle is GTG on codon 7; standard hemoglobin is GAG on codon 7.
dna_seq <- as.character(most_diff_dna_trimmed)
start <- regexpr("ATGGTGC",dna_seq)[1]
codon7 <- substr(dna_seq,start + 18 ,start + 20)
codon7
#This indicates that the individual does carry the sickle mutation.
