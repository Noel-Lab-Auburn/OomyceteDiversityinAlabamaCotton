#!/bin/bash -login

############################################
#     Bioinformatic pipeline               #
#                   #
# Running NBC algorithm through DADA2 with R    #
############################################ 

module load R

R CMD BATCH dada2_assigntax_NBC.R