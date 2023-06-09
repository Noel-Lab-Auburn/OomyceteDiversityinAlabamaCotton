[![DOI](https://zenodo.org/badge/454655246.svg)](https://zenodo.org/badge/latestdoi/454655246)

This repository is under a [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) license. This means it is free to use with citation. Please cite the zenodo DOI above. Full license information can be found in the [LICENSE](LICENSE) file


# Oomycete Diversity in Alabama

The objective of this project was to describe patterns in diversity of oomycetes associated with cotton in Alabama. Oomycetes are major pre- and post-emergent damping off pathogens and no extensive surveys of the species present has been done in Alabama.

We used poisson models to model richness data across a soil texture gradient and permanovas to associate changes in composition with soil textures. 

Then given the results of the survey we did seed rot analysis to determine which species were capable of causing seed rot in the lab. 

Description of R scripts and links:

The two R scripts that generate figures and analysis for this paper are the [Pathogenicity_2.R](Pathogenicity_2.R) file and the [BiogeographyModeling.R](BiogeographyModeling.R) file. 

For the biogeography modeling portion, a phyloseq object was used, and this can be generated from scratch using the species count (OTU) table, Taxonomy Table, and metadata file, or can be directly loaded into R using the oomycete_phyloseq_final.RDS file. The 2023-11-18_IsolateCollection.csv file is all of the isolates in long format and was used to create the species count table. 



