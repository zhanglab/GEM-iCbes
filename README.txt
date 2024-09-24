# This is the git repo for the Genome-scale metabolic model of
Cladicellulosiruptor bescii (GEM-iCbes)

Details related to the metabolic simulations can be found follow the reference:
Zhang K, Zhao W, Rodionov DA, Rubinstein GM, Nguyen DN, Tanwee TNN, Crosby J,
Bing RG, Kelly RM, Adams MWW, Zhang Y. 2021. Genome-scale metabolic model of
Caldicellulosiruptor bescii reveals optimal metabolic engineering strategies
for bio-based chemical production.
mSystems 6:e01351-20. https://doi.org/10.1128/mSystems.01351-20.

This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License.](https://creativecommons.org/licenses/by-nc-nd/4.0/)

# This repository includes the following files:

## README.txt: this file, list of files and directories in the repo

## readme.md: a markdown tutorial on the model simulations

## Genome information
  # Cbes_genome/: directory contains the genome annotation and NCBI protein ID
									mapping of C. bescii strain DSMZ 6725

## Model related
	# model.yaml: model definition YAML file
	# compounds.yaml: compounds database
	# reactions.yaml: reactions database for the wild-type strain (DSMZ 6725)
	# WT_model_def.tsv: reaction list for the wild-type strain (DSMZ 6725)
	# engineered_reactions.yaml: reactions in the engineered strain (MACB1062)
	# EX_DG25.tsv: default exchange constraints of the medium DG25
	# EX_modified_DSMZ516_v1: default exchange constraints of the
														modified DSMZ 516 medium v1
	# EX_modified_DSMZ516_v2: default exchange constraints of the
														modified DSMZ 516 medium v2

## Model simulation scripts
  # scripts/: python and R scripts for the simulation of WT and engineered
							models and the visualization of model results


## additional_files
  # example files used in the simulations described in readme.md
