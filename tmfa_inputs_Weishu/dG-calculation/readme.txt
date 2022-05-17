# This folder is for the standard delta-G calculation for the reactions, based on the cpe structures in Cbes model. 

###############################################################

# (1) group-contribution 

#########

# (1-1) environment set up

## group-contribution: It's a fork in my git repo from "https://github.com/bdu91/group-contribution".

git clone https://github.com/weishuzhao/group-contribution.git


## need a environment of python 2.7 and other softwares:

conda create python=2.7 -n group-contribution
Source activate group-contribution

	# pip install numpy>=1.13.3 scipy>=1.0.0 scikit-learn>=0.19.1 pandas>=0.21.0
	# conda install -c conda-forge openbabel=2
	# conda install -c rdkit rdkit
	# pip install -U lmfit sympy
	# export PATH=$PATH:"/Applications/MarvinSuite/bin/"


# (1-2) test the method by using the code test_group_contribution.py in the group-contribution folder:

cd group-contribution

python test_group_contribution.py
	# Three case studies were provided:
		## [case 1] the same as Keith used to test; 
		## [case 2] R00428, two compounds have gc_db mapping, 1 compound (C05922) used different formats of SMILES (from ModelSeed, converted canonical SMILES, canonical SMILES of the un-charged version from PubChem, and Isomeric SMILES of the uncharged version from PubChem); 
		## [case 3] ATPSYN, to test if any changes happened if all four compounds from gc_db, or some of them using SMILES.



# (1-3) prepre gc_db mapping and SMILES information for compounds

source activate group-contribution

mkdir Cbes-gc-inputs

## Detailed data are shown in the "01-cpd_structure" sheet of table: https://docs.google.com/spreadsheets/d/1Cn9HcrJahVDtvcYK_gCl-CEk-s1azxKHvmtSHOPLBmw/edit?usp=sharing 

cd Cbes-gc-inputs

nano Cbes-gcdb.tsv
nano Cbes-smiles.tsv

cd .. 


# (1-4) calculate the standard delta-G for the Cbes mode (with engineered ethanol produced reactions AdhE and Rnf)

mkdir Cbes-gc-results

## model in the branch "clean-model-ethanol"

# pH 7, IS 0

python test-gc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --gc-db Cbes-gc-inputs/Cbes-gcdb.tsv --smiles Cbes-gc-inputs/Cbes-smiles.tsv --pH 7.0 --IS 0 > Cbes-gc-results/Cbes-dG_pH7_IS0.tsv

# pH 7, IS 0.25

python test-gc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --gc-db Cbes-gc-inputs/Cbes-gcdb.tsv --smiles Cbes-gb-inputs/Cbes-smiles.tsv --pH 7.0 --IS 0.25 > Cbes-gc-results/Cbes-dG_pH7_IS0.25.tsv


# pH 7, IS 0.036

python test-gc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --gc-db Cbes-gc-inputs/Cbes-gcdb.tsv --smiles Cbes-gc-inputs/Cbes-smiles.tsv --pH 7.0 --IS 0.036 > Cbes-gc-results/Cbes-dG_pH7_IS0.036.tsv


# pH 4, IS 0

python test-gc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --gc-db Cbes-gc-inputs/Cbes-gcdb.tsv --smiles Cbes-gc-inputs/Cbes-smiles.tsv --pH 4.0 --IS 0 > Cbes-gc-results/Cbes-dG_pH4_IS0.tsv


# pH 4, IS 0.036

python test-gc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --gc-db Cbes-gc-inputs/Cbes-gcdb.tsv --smiles Cbes-gc-inputs/Cbes-smiles.tsv --pH 4.0 --IS 0 > Cbes-gc-results/Cbes-dG_pH4_IS0.036.tsv



cd ..



#################################################################


# (2) component-contribution

## equilibrator-api is a fork repo: https://gitlab.com/weishuzhao/equilibrator-api.git (original version from: https://gitlab.com/equilibrator/equilibrator-api.git)



source activate component-contribution


mkdir results

## pH 7, IS 0

python ../../../Script_Cbes/TMFA-data_preparsion/script/test-cc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --outfile results/cbes-cc-dG_pH7_IS0.tsv --i 0 --ph 7.0

## pH 7, IS 0.25M (reported for standard aqueous conditions, default set up in component-contribution method)

python ../../../Script_Cbes/TMFA-data_preparsion/script/test-cc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --outfile results/cbes-cc-dG_pH7_IS0.25.tsv --i 0.25 --ph 7.0


## pH 7, IS 0.036M (calculated based on Cbes medium modified DSMZ516 v2)

python ../../../Script_Cbes/TMFA-data_preparsion/script/test-cc-predict-cbes.py --model ../../../GEM_iCbes/GEM-iCbes/model.yaml --outfile results/cbes-cc-dG_pH7_IS0.036.tsv --i 0.036 --ph 7




# combined the results of pH 7 and IS 0.036 M by using both group-contribution and component contribution:

python ../../../Code\ Notes/Script/multicomb.py -l rxn_list.tsv -r dG-gc-pH7_IS0.036.tsv --lci 1 --rci 1 -l rxn_list.tsv -r dG-cc-pH7_IS0.036.tsv --lci 1 --rci 1 -o rxn_list-combined_gc_cc_pH7_IS0.036.tsv





