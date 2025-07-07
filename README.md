# Repository - Coral Microbiome Ontogeny

This repository contains the scripts used to analyze data and create figures for the manuscript "Conspicuous coral-associated bacteria are vertically transmitted from parental gametes"

Raw sequencing data are deposited in the NCBI Sequence Read Archive (SRA) under BioProject [PRJNA1113312](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1113312)


-The script `Larvae_QC.R` was used to identify and remove contaminant ASVs as well as to remove any ASV assigned to chloroplast or mitochondria or ASVs with only zeros. The script also renames sample names and exports a new ASV table.

-The script `Larvae_ANCOMB_barplots.R` was used to plot the most abundant bacterial groups across samples

-The script `Larvae_alpha_diversity.R` was used to calculate alpha diversity and run statistical comparisons using linear models 

-The script `Larvae_beta_diversity.R` was used to calculate beta diversity and run statistical comparisons using PERMANOVAs

-The script  `Larvae_ANCOMB.R` was used to run the differentially abundandance analysis 

-The script  `Larvae_ANCOMB_barplots.R` was used to plot ANCOM-BC results

-The script `Larvae_UpSet_plots.R` was used to plot ASVs that were overlapping across samples

-The script `Larvae_alluvial.R` was used to plot the relative contribution of differentially abundant ASVs across samples
