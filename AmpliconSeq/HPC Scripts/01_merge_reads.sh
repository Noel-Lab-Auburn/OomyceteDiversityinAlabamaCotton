#!/bin/bash
#


##################################
#     Bioinformatic pipeline     #
#     ITS oomycete amplicon      #
#                                                                #
#     ----------------------     #
#        merging reads       #
#                                                                #
##################################

# This script is run from the scripts directory in the scratch. It loops over the names in samples.txt. samples.txt should be in the scripts directory


mkdir merged

#  load the module
module load vsearch

for sample in $(cat samples.txt)
do

    echo "On sample: $sample"
    
    vsearch -fastq_mergepairs ~/noel_shared/Demultiplexed_2023_Oomycetes_NCST2023_AL_Soils/${sample}_R1_001.fastq.gz -reverse ~/noel_shared/Demultiplexed_2023_Oomycetes_NCST2023_AL_Soils/${sample}_R2_001.fastq.gz -fastqout merged/${sample}_merged.fastq -fastq_maxdiffs 20 
    
done

