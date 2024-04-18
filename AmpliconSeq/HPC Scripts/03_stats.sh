#!/bin/bash
#


##################################
#     Bioinformatic pipeline     #
#     ITS oomycete amplicon sequences     #
#                                                                #
#     ----------------------     #
#        stats       #
#                                                                #
##################################

module load vsearch
mkdir stats

cat trimmed/*.fastq > trimmed.fastq

vsearch -fastq_stats trimmed.fastq -log stats/stats_results.txt
