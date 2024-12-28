# PhyloBiomeTurnover
Repository associated with: "Phylogenetic turnover in biomes is determined by environmental harshness, geographic distance, and climate stability"

*Manuscript Authors: Cecile Renfro, Lydia Morley, Daniel Spalink*

*Contact: crenfro\@tamu.edu*

This archive contains the code required to replicate the analyses presented in the main article. The analyses presented in the main article use a large occurrence dataset previously curated and archived by coauthor Lydia Morley. The code in this repository is subsequent to Lydia Morley's code and documentation for occurrence data downloading and manipulation. At the following link, you will find information on the downloading and cleaning of the large dataset, as well as a sample dataset for local machines. Please follow scripts 1, 2A, 2B, and 3 at the following link before running code associated with this repository: <https://github.com/lydMor/MassiveOccurrenceDataManipulation>. Note script 4 is for community matrix gridding which, for biomes, should be done following the script corresponding to this repository.

This repository includes resources for our manuscript pertaining to data manipulation and analyses within and among biomes. It includes an R script ("PhyloBiomeTurnover_20241218.R") that details the data manipulation and analysis of vascular plant biodiversity within and among biomes to understand how the conditions of a biome shape its assembly processes. This R script includes all methods subsequent to the MassiveOccurrenceDataManipulation repository, and the cleaned sample dataset is the input dataset in this script. This script includes assigning occurrences to biomes, biome-specific grid cell delineation, calculating and plotting diversity metrics, GDMs, and more.

The dependencies required to run this script are located at the beginning of the script. Note that while the code is able to run from start to finish on a local machine using a small dataset, it is recommended to break the code into manageable chunks when using a large dataset. The script outputs many files throughout to simplify this process.

This repository also includes two additional R scripts ("phylogdm_functions_sim.R" and "taxgdm_functions_sim.R"), which correspond to the modified phylogenetic and taxonomic formatsitepair functions, respectively. These scripts set up the functions to create a site-pair table for input into the gdm() function in the *gdm* package, but are modified to calculate phylogenetic or taxonomic Simpson's beta diversity. 

Raster files associated with predictor variables (except 19 Worldclim variables, linked below) are included in this repository for simplicity. Please see manuscript text for citations.

The analyses with the full GBIF dataset were run on a supercomputer using R version 4.2.2 (2022-10-31) "Innocent and Trusting".

The GBOTB.tre Smith and Brown (2018) seed plant phylogeny can be accessed at this repository: <https://github.com/FePhyFoFum/big_seed_plant_trees/releases>

The Nita et al. (2022) fern tree phylogeny v1.5.0 can be accessed at this repository: <https://fernphy.github.io/downloads.html>

The Dinerstein et al. (2017) biome shapefile can be accessed here: <https://ecoregions.appspot.com/> \> About \> Download data \> Shapefile (150mb zip)

Worldclim bioclimatic raster brick can be accessed here: https://www.worldclim.org/data/worldclim21.html > Bioclimatic variables > bio 30s
