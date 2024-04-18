#!/bin/bash
#


##################################
#     Bioinformatic pipeline     #
#     Oomycete ITS     #
#                                                                #
#     ----------------------     #
#        stripping primers       #
#                                                                #
##################################


mkdir trimmed

#  load the module
module load anaconda/3-2020.02

for sample in $(cat samples.txt)
do

    echo "On sample: $sample"
    
    cutadapt -g GAAGGTGAAGTCGTAACAAGG -a AGCGTTCTTCATCGATGTGC -f fastq -n 2 -m 20 --discard-untrimmed --match-read-wildcards merged/${sample}_merged.fastq > trimmed/${sample}_trimmed.fastq

done