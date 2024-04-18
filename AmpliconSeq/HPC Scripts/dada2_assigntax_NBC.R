library(dada2)
library(Biostrings)

taxonomy.file.path <- "~/noel_shared/db_fungi/sh_general_release_dynamic_s_all_25.07.2023_mockseqadded.fasta"
otus.file.path <- "clustered/otus.fasta"

# Fasta
FASTA.otus <- readDNAStringSet(otus.file.path, format="fasta", seek.first.rec=TRUE, use.names=TRUE)

taxa <- assignTaxonomy(FASTA.otus, taxonomy.file.path, multithread=TRUE, tryRC = TRUE)

saveRDS(taxa, file = "taxa_out_DADA2_NBC.rds")
