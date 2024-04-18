#!/bin/bash 

############################################
#     Bioinformatic pipeline                #
#     ITS Oomycete sequenicng             #
#DEREPLICATION, CLUSTERING, CHIMERA REMOVAL#
############################################

#HAVE TO INSTALL USEARCH BEFORE USING

module load vsearch
mkdir clustered

# dereplication 
vsearch --derep_fulllength filtered/filtered_trimmed.fasta --output filtered/uniques.fasta -sizeout

# de-noising (error correction), output is zero-radius OTUs
#usearch -unoise3 output/clustered/uniques_R1.fasta -tabbedout output/clustered/unoise_zotus_R1.txt -zotus output/clustered/zotus_R1.fasta

# clusters OTUs based on traditional 97% identity 
usearch -cluster_otus filtered/uniques.fasta -minsize 2 -otus clustered/otus.fasta -uparseout clustered/uparse_otus.txt -relabel OOTU_ --threads 20

# useful links
#http://www.drive5.com/usearch/manual/unoise_pipeline.html
#http://www.drive5.com/usearch/manual/faq_uparse_or_unoise.html
#http://www.drive5.com/usearch/manual/cmd_otutab.html
#http://www.drive5.com/usearch/manual/upp_labels_sample.html
#http://drive5.com/usearch/manual/bugs.html
#http://drive5.com/usearch/manual/support.html
