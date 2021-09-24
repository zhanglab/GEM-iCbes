# iCbes model, Genome-scale metabolic model of *Cladicellulosiruptor bescii*  

Files in this repository represent the genome-scale model of *Cladicellulosiruptor bescii* (*C.bescii*). The model is first reported in the following publication:

* Zhang K, Zhao W, Rodionov DA, Rubinstein GM, Nguyen DN, Tanwee TNN, Crosby J, Bing RG, Kelly RM, Adams MWW, Zhang Y. 2021. Genome-scale metabolic model of *Caldicellulosiruptor bescii* reveals optimal metabolic engineering strategies for bio-based chemical production. *mSystems* 6:e01351-20. https://doi.org/10.1128/mSystems.01351-20.

## List of files  
* model.yaml:
  > model definition YAML file  
* compounds.yaml:
  > definition of metabolites
* reactions.yaml:
  > definition of reactions in the model based on genome annotations of the wild-type *C. bescii* strain DSM 6725  
* WT_model_def.tsv:
  > list of reactions included in the wild-type *C. bescii* model  
* engineered_reactions.yaml:
  > definition of reactions related to the simulation of engineering designs in C. bescii
* EX_DG25.tsv:
  > default exchange constraints of the DG25 medium  
* EX_modified_DSMZ516_v1:
  > default exchange constraints of the modified DSMZ 516 medium in version 1  
* EX_modified_DSMZ516_v2:
  > default exchange constraints of the modified DSMZ 516 medium in version 2  
* scripts/:
  > directory containing python scripts used for model simulations
* additional_files/:
  > directory containing additional input files that are used for model simulations.


# Model simulations
Model simulations were performed using the open source software [PSAMM](https://zhanglab.github.io/psamm/) [1]. Additional python & R scripts were used to set up computational growth simulations with the *C. bescii* GEM. A detailed protocol of model simulations is provided below.

## Materials  
* **PSAMM** is an open source software that is designed for metabolic model curation and analysis. To install `PSAMM` and the associated requirements, you can reference the [Installation and Materials](https://psamm.readthedocs.io/en/latest/tutorial/psamm-install.html) section of PSAMM documentation.
* **R needs to be installed**. If you haven't install R, following the instruction [here](https://www.r-project.org/) to install it.
* You should see python scripts in a folder called **scripts** in this repository.
* A folder called **additional_files** contains some additional input files that will be used to run the commands in this
tutorial.

Create a folder named **simulations** and navigate to that folder, you will run simulation commands in this directory and  save the simulation outputs here:  
```shell
(psamm-env) $ mkdir simulations
(psamm-env) $ cd simulations
```

## Optimizing yields of growth or products   
`scripts/fva_exchange_constrain.py` and `scripts/FBA_MOMA_exchange_constrain.py` are used for model validation. Both the growth yields and the product generation were compared to experimental measurements. Specific constraints (e.g. exchange flux bounds, reaction flux bounds, gene knockouts) can be applied via argument inputs of the scripts. More information on parameter options can be found in the help message using the following commands:  
```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py -h
(psamm-env) $ python ../scripts/FBA_MOMA_exchange_constrain.py -h
```

Here, we use `fva_exchange_constrain.py` as example to show the detailed usage of the scripts. **Check the model.yaml and make sure you use correct medium file when running simulations below**

### **Basic usage**
To run `fva_exchange_constrain.py`, a `--model` option should be specified to provide a path to the `model.yaml` file:
```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py --model ../model.yaml
```
By default, this command will apply flux variability analysis (FVA), which identifies the optimized biomass yield and the simulate the maximum and minimum flux of each reaction when the biomass yield is reaching that optimal. The input settings are specified in the `model.yaml` file. The output contains five columns, indicating the following information: *reaction id*, *minimum flux*, *maximum flux*, *flux range*, *reaction equation*. An example is shown below:
```
R01210  0.0     0.0     [0.0, 0.0]      |3-Methyl-2-oxobutanoic acid[c]| + CoA[c] + NAD+[c] => 2-Methylpropanoyl-CoA[c] + CO2[c] + NADH[c]
Ram106B 0.0     0.0     [0.0, 0.0]      |Rhamnose oligosaccharides(n=2)[c]| + H2O[c] => (2) L-Rhamnose[c]
sink_biomass    0.0     0.0     [0.0, 0.0]      Biomass[c] =>
...
```
In the above example, the optimized biomass is 0. This is because the default model provides no carbon source in the medium.

### **Specify the simulation approach**
Specific simulation approaches can be applied by setting the `--method` option, which can be assigned to one of the followng: `fba`, `fva` (for running `fva_exchange_constrain.py`), or `fba`, `lin_moma`, `lin_moma2`, `moma`, `moma2` (for running `FBA_MOMA_exchange_constrain.py`). You can find more information about these methods in the [PSAMM documentation](https://psamm.readthedocs.io/en/latest/commands.html). An example is shown below:

```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py --model ../model.yaml --method fba
```
The runs the `fba` method with `fva_exchange_constrain.py`. It will generate the following output, with the *reaction id*, *reaction flux*, *reaction equation* specified in the corresponding columns:
```
R01210  0.0     |3-Methyl-2-oxobutanoic acid[c]| + CoA[c] + NAD+[c] => 2-Methylpropanoyl-CoA[c] + CO2[c] + NADH[c]
Ram106B 0.0     |Rhamnose oligosaccharides(n=2)[c]| + H2O[c] => (2) L-Rhamnose[c]
sink_biomass    0.0     Biomass[c] =>
...
```

### **Add exchange constraints**
`--exchange-list` is the argument used to specify exchange constraints (i.e. media uptake and/or product diffusion) to the model. This is done by providing a TSV-file with three columns that specify: *compound id*, *lower bound*, *upper bound*. Compounds in the TSV-file are assigned to the `extracellular_compartment` by default. An example is given in below:

Exchange list file:
```python
# lactate (C00186) exchange with production limit of at least 22.004
C00186  22.004  1000
# glucose (C00031) exchange with uptake limit of -16.83
C00031  -16.83  1000
...
```

Command-line:
```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py --model ../model.yaml --exchange-list ../additional_files/B0-G-b.tsv --method fba
```    

### **Specify objective reaction(s) and add metabolic reaction constraints**
The argument `--objective` is used to specify the objective function. If this option is not assigned, the default biomass function designated in the `model.yaml` file will be used as objective for simulations. The argument expect to receive a comma-separated list of reaction IDs in the model. If only one reaction is listed, then that reaction will be optimized; if more than one reaction is listed, then a linear combination of reaction fluxes across all reactions are optimized. A related  argument `--rule` is used to constrain the ratio of reaction fluxes when multiple reactions was listed in the `--objective` option. The `--rule` option is defined by a string like "**ratio1,reaction1_ID,ratio2,reaction2_ID**", which means **flux(reaction1) : flux(reaction2) = float(ratio1) : float(ratio2)**. Hence, it maintains a constant flux ratio between two reaction fluxes. Another relevant argument is the `--rxn-constrain` option, which is followed by a TSV-file containing three columns: reaction id, lower bound (of reaction flux), upper bound (of reaction flux). This can be used to constrain the biomass yield to a non-zero value for the optimization of product yields. For example:  

```python
# reaction sink_biomass (for growth yields) with a minimum biomass yield of 0.5251
sink_biomass    0.5251  1000
```  

Below is an example command on the optimization of product yields with specific contraints on the product ratios:  
```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py --model ../model.yaml --exchange-list ../additional_files/B0-G-p.tsv --objective EX_C00186[e],EX_C00033[e],EX_C00022[e] --rule 1,EX_C00186[e],0.33133,EX_C00033[e] --rule 1,EX_C00186[e],0.006,EX_C00022[e] --rxn-constrain ../additional_files/B0-G-p_biomass_constraint.tsv --method fba
```
In the above command, the objective function is sum of lactate (EX_C00186[e]), acetate(EX_C00033[e]) and pyruvate (EX_C00022[e]) production, with the constraints of acetate_production : lactate_production = 0.33133 and pyruvate_production : lactate_production = 0.006, and growth yields (flux of sink_biomass) equal to or larger than 0.5251 g/L.

### **Simulating the engineered strains**
The default model represents wild-type *C. bescii*, deleting reactions associated with specific genes or adding engineered reactions is needed in order to simulate mutant  or engineered strains. The arguments `--gene` and `--addrxn` are used for this purpose. `--gene` accepts a space-delimited list of gene IDs, and remove all reactions that require any of the genes on the list. `--addrxn` is followed by a space-delimited list of reaction IDs (**Note:** these reactions should be in the reactions.yaml but not in the file WT_model_def.tsv). An example command for using the `--gene` and `--addrxn` options is given below:  
```shell
(psamm-env) $ python ../scripts/fva_exchange_constrain.py --model ../model.yaml --gene Athe_1382 Athe_1918 --addrxn R00228 R00754 Rnf_Na --objective EX_C00469[e] --exchange-list ../additional_files/E1-C-p.tsv --rxn-constrain ../additional_files/E1-C-p_biomass_constraint.tsv --method fba
```
The command above simulate the engineering of ethanol-producing *C. bescii*. It deletes two reactions that are catalyzed by enzymes encoded by the genes Athe_1382 (R08231, R01870) or Athe_1918 (R00342, R00703), and add three engineered reactions (R00228, R00754, Rnf_Na) from the *C. bescii* model. Additional constraints to the compound exchange is defined in `E1-C-p.tsv`, and constraint to the biomass flux is defined in `E1-C-p_biomass_constraint.tsv`. The ethanol production (EX_C00469[e]) was optimized using the `fba` method.

## Carbon utilization test
`fva_change_CarbonInput.py` was used to simulate the growth of *C. bescii* on different sole carbon sources, multiple carbon sources can be tested simultaneously using this script. It requires a TSV input file with four columns: *compound id*, *lower bound of exchange constraints*, *upper bound of exchange constraints*, *compound name*. An example of the input file is given below:

```python
# L-arabinose with uptake limit of 1
C00259  -1      1000    L-arabinose
C00185  -1      1000    cellobiose
C00760_kz       -1      1000    crystalline cellulose
C00095  -1      1000    fructose
C00124  -1      1000    galactose
C00031  -1      1000    glucose
...
```
An example command:  
```python
python ../scripts/fva_change_CarbonInput.py --model ../model.yaml --objective sink_biomass --exchange-list ../additional_files/carbon_sources_1mM.tsv
```
The output contains nine columns, as shown in the example output below:
```
carbon_source   compound_name   formula objective       lower_bound_of_objective        upper_bound_of_objective        lower_bound_of_carbon upper_bound_of_carbon carbon_used_up?
C00259  L-Arabinose     C5H10O5 sink_biomass    0.064168        0.064168        -1.0    -0.999990834175 Yes
C00185  Cellobiose      C12H22O11       sink_biomass    0.219966        0.219966        -1.0    -0.999995651544 Yes
C00760_kz       Cellulose(n=218)        H2O(C6H10O5)218 sink_biomass    0.912923        0.912923        -1.0    -0.038084389863 No
C00095  D-Fructose      C6H12O6 sink_biomass    0.110003        0.110003        -1.0    -0.999997327495 Yes
C00124  D-Galactose     C6H12O6 sink_biomass    0.082502        0.082502        -1.0    -0.999994297259 Yes
...
```

## Minimal network analysis
The minimal network (MN) analysis was used for the identification of essential functions for the optimization of bioproduct production in *C. bescii*, such as ethanol. The protocol in applying the MN analysis is (using ethanol as an example):
* **(1)** Perform three sets of 1,000 randomsparse simulations with ethanol production constrained to 0% (no-ethanol), 50% (half-maximum), and 99.99% (maximum) of the predicted maximum, respectively.
* **(2)** Classify all reactions into three categories: **(i) core-essential** set, representing reactions required in all minimal networks of a simulation condition; **(ii) conditionally-essential** set, representing reactions required in some but not all minimal networks; **(iii) non-essential** set, representing reactions not used in any minimal networks.
* **(3)** Compare results at different ethanol threshold to identify the essential functions for each condition.
* **(4)** Perform a convergence test for each MN simulation to confirm that the classification of core-essential, conditionally-essential, and non-essential reactions were stabilized after 1,000 random simulations.

### **E1M strain specification**
The MN analysis was performed with the **E1M** strain, with a genotype of *adhE<sup>+</sup>rnf_Na<sup>+</sup>*. This requires a modification on the default (WT) model to make it represent E1M.  
* (1) E1M strain setup: modify the `model.yaml` file to add the engineered reactions (by commenting out the `WT_model_def.tsv` specification assigned under the `model:` option), keep only the `EX_modified_DSMZ516_v2.tsv` as the exchange settings, and add a `limits.yaml` file under the `limits:` option to specify the biomass constraints. The resulting `model.yaml` file is shown in below:

    ```python
    name: iCbes, Genome-scale metabolic model of Cladicellulosiruptor bescii
    biomass: Biomass_Cbescii
    default_flux_limit: 1000
    default_compartment: c
    extracellular: e
    compartments:
    - id: c
    adjacent_to: e
    name: Cytoplasm
    - id: e
    adjacent_to: c
    name: Extracellular
    compounds:
    - include: compounds.yaml
    reactions:
    - include: reactions.yaml
    - include: engineered_reactions.yaml  # added to simulate the C. bescii engineered ethanol strain MACB1062
    exchange:
    #- include: EX_DG25.tsv # used to test the growth and product yields of wild-type and dLdh strain grown on glucose or fructose, as well as the carbon utilization for model validation
    #  format: tsv
    #- include: EX_modified_DSMZ516_v1.tsv # used to test the growth yields of wild-type strain on cellulose
    #  format: tsv
    - include: EX_modified_DSMZ516_v2.tsv # used to test the ethanol production of the engineered ethano strain MACB1062 on cellulose
      format: tsv
    #model:
    #- include: WT_model_def.tsv
    limits:
    - include: limits.yaml
    ```

* (2) Set cellulose as the sole carbon source in `EX_modified_DSMZ516_v2.tsv`:
    ```python
    ...
    C00760_kz	e       -0.3337  1000
    ...
    ```

* (3) Constrain sink_biomass in `limits.yaml` based on experimental data:
    ```python
    - reaction: sink_biomass
      lower: 0.41
      upper: 1000
    ```

After the above model editing, change directory to the `simulations/` folder. Then we should be ready for making the MN simulations:  
```shell
(psamm-env) $ cd simulations
```

### **Reaction randomsparse**
This randomsparse simulation performs successive random deletions of metabolic genes or reactions while maintaining a user-defined minimum flux on an objective function. An example command for randomsparse simulation with reaction deletions is given below:
```shell
(psamm-env) $ psamm-model --model <path to model.yaml> randomsparse --objective <objective reaction ID> --type reactions <threshold>
```
* **Note:** In our examples below, the last parameter for the command, \<threshold\>, is given as a percentage value that will be applied to the maximum objective flux. The successive deletion will stop when no further deletion can be made while still maintaining the objective flux to above or equal to the threshold.  

For example, the following command runs reaction randomsparse with an optimized ethanol production of at least 99.99% of the maximum:  
```shell
(psamm-env) $ psamm-model --model ../model.yaml randomsparse --objective EX_C00469[e] --type reactions 99.99%
```

The output contains two columns: *reaction id*, composed of values 0 or 1, with 1 indicating an essential reaction in a random minimal network and 0 indicating a deleted reaction in the minimal network. An example output can be seen as follows:
```
INFO: Essential reactions: 329/759
INFO: Deleted reactions: 430/759
12DAG3PS	1
12PD_TP	0
14GLACANASE	0
16GLACANASE	0
ACt6	1
ANTIMt1	0
ARSNAt1	0
ARSt1	0
ATPSYN	1
AXE	0
AbfA	0
...
```

For more information of randomsparse, see [Random Minimal Network Analysis](https://psamm.readthedocs.io/en/latest/tutorial/constraint_based.html#random-minimal-network-analysis) in the PSAMM documentation.

### **Convergence test**
A convergence test can be performed over a high number (e.g. 1,000) of randomsparse simulations with the same model definition. The simulation results should be combined into a single TSV-file (For 1,000 random simulations, this file should contain 1001 columns: one for the reactions IDs and 1,000 collumns of corresponding 0/1 values for the  randomsparse simulation results). An example of the combined TSV file is in `../additional_files/E1M-mn_99.99_combine_all_repeats.tsv `.

`sufficiency_test.R` is used to perform the convergence test. The general usage is:
```shell
(psamm-env) $ Rscript ../scripts/sufficiency_test.R <path to combined TSV result file> 0 1000 <label of the objective being optimized> <Number of Genes/Reactions>
```

An example command is shown below:  
```shell
(psamm-env) $ Rscript ../scripts/sufficiency_test.R ../additional_files/E1M-mn_99.99_combine_all_repeats.tsv 0 1000 "Ethanol production" "Number of reactions"
```
This command will create two files as outputs:
* A plot for the visualization of convergence.
* A TSV-file indicating the number of core-essential, conditionally-essential and non-essential reactions with the increasing number of random simulations.

### **Classify reactions**
`randomsparse_repeats_classfication.py` is used to classify the reactions into three categories: core-essential, conditionally-essential, and non-essential. An example command is in below:
```shell
(psamm-env) $ python ../scripts/randomsparse_repeats_classfication.py --input ../additional_files/E1M-mn_99.99_combine_all_repeats.tsv --repeat 1000 --out1 E1M-mn_99.99_combine_all_repeats_classification.txt
```

The output, `E1M-mn_99.99_combine_all_repeats_classification.txt` contains four columns: *reaction_id*, *classification*, *number of simulations where this reaction is essential in the minimal network*, *number of simulations where this reaction is not essential in the minimal network*. An example is shown below:
```
reaction_id classification  in_minimal_network      not_in_minimal_network
R09726  Non-essential   0       1000
R10244  Non-essential   0       1000
R11398_kz       Non-essential   0       1000
R04591  Core    1000    0
R04945  Flexible        353     647
R04944  Non-essential   0       1000
R04946  Flexible        353     647
...
```

## Gene knockin-knockout simulations
Gene knockin-knockout simulations were used to evaluate the effect of specific engineering plans on the optimized ethanol production in *C. bescii*. `fva_gene_knockout_knockin.py` was designed to perform this kind of simulation. Two inputs are required:

* **A Knockout_knockin table:** A 3-column TSV-file that define the output file name of a simulation, the knockout genes, and/or engineered knockin reactions. Each line in the file represent a single simulation. For example:
    ```python
    # no gene knockout and no reaction addition, name of simulation result is "base-strain_FVA.tsv"
    base-strain_FVA.tsv     --      --
    # knock out genes encoding bifurcating hydrogenase (BF-H2ase), no reaction addition
    dBFH2ase_FVA.tsv        Athe_1295 Athe_1296 Athe_1297 Athe_1298 Athe_1299       --
    # knock out genes encodig PyrE and lactate dehydrogenase (LDH), add reaction Mrp
    Mrp+dPyrE+dLDH_FVA.tsv  Athe_1382 Athe_1918     Mrp
    ```
* **An output file path:** path to an output table summarizing the simulation conditions. Columns of the output  table indicates: *flux of the objective function*, *flux range of the objective function*, the exchange flux ranges of *acetate*, *pyruvate*, *lactate*, *acetoin*, *CO<sub>2</sub>*, *H<sub>2</sub>*, and the range of biomass flux under each corresponding condition listed in the knockout_knockin table.

Besides, you can define an objective function (the default objective is the biomass reaction defined in model.yaml), the simulation method (`fva` or `fba`, with the default method as `fva`) in command line when running the gene knockin-knockout simulations.

**Before performing the simulation, set up the model according the descriptions in the *E1 strain specification* section, uncomment four reactions (Rnf_H, Mrp, R07181 and R00700) in the  `engineered_reactions.yaml` file**  

An example command can be seen as follows:  
```shell
(psamm-env) $ python ../scripts/fva_gene_knockout_knockin.py --model ../model.yaml --objective "EX_C00469[e]" --knock-table ../additional_files/knockout_knockin.list --method fva --decimal 4 --final-print ./01-opt_EtOH_final_result.tsv
```
This example command will export three FVA simulation results called `base-strain_FVA.tsv`, `dBFH2ase_FVA.tsv` and `Mrp+dPyrE+dLDH_FVA.tsv`, respectively. As well as a file called `01-opt_EtOH_final_result.tsv`, which summarizes fluxes of targeted reactions under the three specified knockout-knockin conditions.

## Robustness-fva analysis
The robustness-fva analysis is used to evaluate how the flux ranges of targeted metabolic reactions would change under increasing fluxes of another reaction. It calculates the maximum and minimum fluxes of a reaction along each step of the varying flux of another reaction, which setting the optimal of an objective function (e.g. biomass).

### **Run robustness-fva**
The tutorial below provides examples on analyzing the fluxes of targeted reactions with varying ethanol production, from zero to the maximum in engineered *C. bescii* models. **The model should be configured according to specified in the *Gene knockin-knockout simulations* section.**

#### **Basic use of robustness-fva**
The basic command can be run as the following:
```shell
(psamm-env) $ python ../scripts/robust-fva.py --model ../model.yaml --vary EX_C00469[e] --objective sink_biomass --steps 100
```
The reaction that is maximized and minimized is designated with the `--objective` option. The reaction that is varied is designated with the `--vary` option. The number of steps in the robustness-fva simulation is assigned with the `--steps` option.  

The output contains four columns: *reaction ID*, *flux of the varying reaction in each step* (in this case the sink_biomass reaction), *the minimum flux of the target reaction in a given step*, and *the maximum flux of the target reaction in a given step*. An example output is given below:  
```
R09726  0.0     0.0     0.0
R10244  0.0     0.0     0.0
R11398_kz       0.0     0.0     0.0
R04591  0.0     0.180787949     0.181086822
R04945  0.0     -1000.0 1000.0
R04944  0.0     0.0     0.0
R04946  0.0     -1000.0 1000.0
R07238  0.0     -1000.0 1000.0
Bgal    0.0     0.0     0.0
...
```
#### **Performing robustness-fva on engineered strains**
Further modifications to the model can be specified with `--gene` and `--addrxn` options when running the robustness-fva simulation. Usage of these two options are the same as described in the *Simulating the engineered strains* section.  

For example, to run the same simulation as above but with the sodium-dependnet Rnf (Rnf_Na) replaced by the proton-dependent Rnf (Rnf_H), one would use the follow command:
```shell
(psamm-env) $ python ../scripts/robust-fva.py --model ../model.yaml --vary EX_C00469[e] --objective sink_biomass --gene engineered_Rnf --addrxn Rnf_H  --steps 100
```
* **Note:** Genes associated with Rnf_Na is labeled as 'engineered_Rnf'; 'Rnf_H' is a reaction ID; Rnf refers to membrane-bound reduced ferredoxin NAD oxidoreductase.

#### **Specify objective and add constraints**
You can customize exchange and metabolic reaction constraints with optional arguments `--exchange-list` and `--rxn-constrain`. And the objective function, specified with option `--objective`, can be the sum of multiple reactions. Furthermore, when there are multiple reactions defined with `--objective`, an optional argument `--rule` can be used to constrain the flux ratio of any two reactions in objective function. For the detailed usage of these four arguments, you can reference the *Add exchange constraints* and *Specify objective reaction(s) and add metabolic reaction constraints* section under the **Optimizing yields of growth or products**.

### **Visualize metabolic fluxes of target reactions**
In order to visualize the flux range of targeted reactions under the varying flux of a reference reaction, the `robust-fva-plot-any-rxn.R` can be used. The general usage of this script is:
```python
Rscript scripts/robust-fva-plot-any-rxn.R <infile> <varylabel> <plotted_reactions> <panel_dimension> <output_name_suffix>
```  
Below are an explanation of each parameter:
* **infile**: the path of robustness-fva result  
* **varylabel**: the x-axis label, based on the varied reaction, e.g. if vary etahnol production, the varylabel could be 'Ethanol production'.  
* **plotted_reactions**: path to a tab-separated two-column input file that specify the reactions to be plotted, the first column lists reaction IDs  (e.g. R01196), the second column lists reaction labels (e.g.POR). An example of the plotted reactions table is shown below:  
    ```
    ATPSYN  ATPSYN
    BF-Nfn  BF-Nfn
    BF-H2ase        BF-H2ase
    MBH     MBH
    EX_C00282[e]    H2
    Rnf_Na  Rnf_Na
    R01061  GAPDH
    R07159  GOR
    ...
    ```
* **panel_dimension**: The number of panels, c(row, col), in the output .png file. For example, '4,2' means there will be four plots per column and two plots per row in the final .png image.  
* **output_name_suffix**: Suffix of the output image file. If the `infile` is set as 'E1M-robustness-fva-result.tsv' and the `output_name_suffix` is set as '-plot', then the final  output image name is 'E1-robustness-fva-result.tsv-plot.png'. By default, the suffix is `-plot-rxn`.    

An example command is shown as the following:
```python
# visualize flux range of targeted reactions
(psamm-env) $ Rscript ../scripts/robust-fva-plot-any-rxn.R ../additional_files/E1M-0-robust.tsv "Ethanol Production (mM)" ../additional_files/rxns4plot-8keyRxns.tsv 4,2 -Cbes-key-enzymes
```
The output will be named "E1M-0-robust.tsv-Cbes-key-enzymes.png", it contains eight plots. In each plot, black (blue) lines with numbers indicate the maximum (minimum) fluxes carried out by specific reactions in the designated model. The area shaded in gray indicates flux ranges of a reaction through the varying steps of the reference reaction.


# References:  
[1] Steffensen JL, Dufault-Thompson K, Zhang Y. PSAMM: A Portable System for the Analysis of Metabolic Models. PLOS  
    Comput Biol. Public Library of Science; 2016;12: e1004732. doi:10.1371/journal.pcbi.1004732.
