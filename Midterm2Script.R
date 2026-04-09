#Dylan Rider
#Midterm 2 script

# installing packages
install.packages("UniprotR")
install.packages("ape")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

# load libraries
library(ape)
library(Biostrings)
library(UniprotR)

# read in the phylogeny tree from iqtree
tree <- read.tree("metazoa_alignment.5k.fasta.treefile")

# root the tree on the porifera branch
rooted_tree <- root(tree, outgroup = c("Plakina_jani", "Grantia_compressa"), resolve.root = TRUE)

# write the rooted tree and plot it
write.tree(rooted_tree, "rooted_metazoa_tree.newick")
pdf("rooted_metazoa_tree.pdf", width = 12, height = 10)
plot(rooted_tree, cex = 0.8)
add.scale.bar()
dev.off()

# split the iqtree support labels into sh-alrt and ultrafast bootstrap
support_raw <- rooted_tree$node.label
support_split <- strsplit(support_raw, "/")
support_table <- data.frame(
  node = seq_along(support_split),
  sh_alrt = suppressWarnings(as.numeric(sapply(support_split, `[`, 1))),
  ufboot = suppressWarnings(as.numeric(sapply(support_split, `[`, 2)))
)

# summarize support values
support_summary <- data.frame(
  metric = c("mean_sh_alrt", "mean_ufboot", "nodes_ufboot_ge_95", "nodes_ufboot_ge_80", "nodes_ufboot_lt_80"),
  value = c(
    mean(support_table$sh_alrt, na.rm = TRUE),
    mean(support_table$ufboot, na.rm = TRUE),
    sum(support_table$ufboot >= 95, na.rm = TRUE),
    sum(support_table$ufboot >= 80, na.rm = TRUE),
    sum(support_table$ufboot < 80, na.rm = TRUE)
  )
)
write.csv(support_table, "tree_support_table.csv", row.names = FALSE)
write.csv(support_summary, "tree_support_summary.csv", row.names = FALSE)

# test monophyly of the major clades from the assignment
lophotrochozoa <- c("Cerebratulus_sp", "Golfingia_vulgaris", "Leptochiton_rugatus", "Macrostomum_lignano", "Novocrania_anomala")
ecdysozoa <- c("Parhyale_sp.", "Priapulus_sp", "Strigamia_maritima", "Tribolium_castaneum")
xenambulacraria <- c("Ptychodera_flava", "Rhabdopleura_sp", "Strongylocentrotus_purpuratus", "Xenoturbella_bocki")
chordata <- c("Homo_sapiens", "Latimeria_chalumnae", "Lepisosteus_oculatus", "Petromyzon_marinus", "Xenopus_tropicalis")
cnidaria <- c("Antipathes_caribbeana", "Plumapathes_pennacea")
porifera <- c("Grantia_compressa", "Plakina_jani")

monophyly_results <- data.frame(
  clade = c("Lophotrochozoa", "Ecdysozoa", "Xenambulacraria", "Chordata", "Cnidaria", "Porifera"),
  monophyletic = c(
    is.monophyletic(rooted_tree, lophotrochozoa),
    is.monophyletic(rooted_tree, ecdysozoa),
    is.monophyletic(rooted_tree, xenambulacraria),
    is.monophyletic(rooted_tree, chordata),
    is.monophyletic(rooted_tree, cnidaria),
    is.monophyletic(rooted_tree, porifera)
  )
)
write.csv(monophyly_results, "clade_monophyly_results.csv", row.names = FALSE)

# read in the gene alignment
sequence_dna <- readDNAStringSet("metazoa_alignment.gene.fasta")

# count dashes and Ns in each sample
gap_summary <- data.frame(
  taxon = names(sequence_dna),
  gaps = vcountPattern("-", sequence_dna),
  Ns = vcountPattern("N", sequence_dna, fixed = TRUE),
  alignment_length = width(sequence_dna)
)
gap_summary$non_gap_non_N <- gap_summary$alignment_length - gap_summary$gaps - gap_summary$Ns
write.csv(gap_summary, "gene_gap_summary.csv", row.names = FALSE)

# find the human sequence
human_sequence <- sequence_dna[names(sequence_dna) == "Homo_sapiens"]

# remove dashes and Ns from the human sequence
human_sequence_clean <- gsub("-", "", as.character(human_sequence))
human_sequence_clean <- gsub("N", "", human_sequence_clean)

# trim sequence to a full codon length and translate it
trim_length <- nchar(human_sequence_clean) - (nchar(human_sequence_clean) %% 3)
human_sequence_clean <- substr(human_sequence_clean, 1, trim_length)
human_dna <- DNAStringSet(human_sequence_clean)
names(human_dna) <- "Homo_sapiens_gene"
writeXStringSet(human_dna, "sequence.fasta")

protein_seq <- translate(human_dna)
names(protein_seq) <- "Homo_sapiens_protein"
writeXStringSet(protein_seq, "protein_sequence.fasta")

# add the best uniprot accession after searching the protein sequence
accessions <- c("P54098")
writeLines(accessions, "uniprot_accessions.txt")

# run uniprot functions once the accession has been filled in
go_info <- GetProteinGOInfo(accessions)
PlotGOAll(GOObj = go_info, Top = 10, directorypath = getwd(), width = 8, height = 5)

# write out additional uniprot tables
names_taxa <- GetNamesTaxa(accessions)
write.csv(names_taxa, "sequence_names_taxa.csv", row.names = FALSE)

