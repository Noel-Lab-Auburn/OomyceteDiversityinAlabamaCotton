#!/bin/bash
#


##################################
#     Bioinformatic pipeline     #
#     ITS oomycete     #
#                                                                #
#     ----------------------     #
#        filtering       #
#                                                                #
##################################


module load vsearch
mkdir filtered

vsearch -fastq_filter trimmed/trimmed.fastq -fastq_maxee 1 -fastq_trunclen 294 -fastq_maxns 0 -fastaout filtered/filtered.fasta -fastqout filtered/filtered.fastq
vsearch -fastq_filter filtered/filtered.fastq -fastq_stripleft 44 -fastaout filtered/filtered_trimmed.fasta -fastqout filtered_trimmed.fastq
