#!/bin/bash 

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences for oomycetes              #
#            MAPPING                       #
############################################

# First use seqtk to convert all the demultiplexed samples into fasta on a loop.
# Use samples.txt like we did for cutadapt

# Load the modules
module load gcc/11.2.0
module load seqtk/1.3-olt7cls
source /apps/profiles/modules_dmc.sh.dyn
module load python/3.3.2
module load anaconda/3-2023.03

mkdir mapping

for sample in $(cat samples.txt)
do

echo "On sample: $sample"
    seqtk seq -a merged/${sample}_merged.fastq > mapping/${sample}_merged.fasta

    # have to replace the beginning of the fasta headers with the file name for mapping. Otherwise we # get one sample with all the read counts, which is not what we want.

    python3 ~/noel_shared/python_scripts/replacefastaheaders_filename.py mapping/${sample}_merged.fasta

done

# have to create one file containing all the reads from the demultiplexed reads
cat *_newheaders.fasta > mapping/demultiplexed_new.fasta


# align the demultiplexed reads back to the now clustered OTUs or ZOTUs (ESV)
module load anaconda/3-2021.11
module load vsearch

vsearch -usearch_global mapping/demultiplexed_new.fasta -db clustered/otus.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_oomycete.txt