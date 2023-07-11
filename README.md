# HOP: the High-throughput Oligogenic Prioritizer

This repository provides the code to run the HOP predictor and reproduce the 
results presented in the paper "Prioritization of Oligogenic Variant Combinations in Exomes"

## 1) Requirements

### Data requirements

The knowledge graph BOCK needs to be downloaded from zenodo (https://doi.org/10.5281/zenodo.7185680) -> BOCK-graphml
extracted and stored in the folder `data/BOCK_data/` under the name `bock.graphml` or the config.ini file should
be updated. 

In order to run the tool, patient exomes need to be annotated with features used by the VarCoPP2.0 predictor. 
Annotated files for the 
1KGP exomes that were used as templates to generate the synthetic exomes, as well as annotated files for the 
OLIDA combinations to be inserted in the 1KGP exomes in order to be able to reproduce the data 
presented in the paper can be downloaded from zenodo : https://doi.org/10.5281/zenodo.8121283
The annotated_data archive should be downloaded and the content extracted in the `data` directory
(`tar -xvf annotated_data.tar.gz`). If you save the data elsewhere, the config.ini file should be updated. 
This archive contains annotated data for 100 exomes (~45GB when extracted).

The annotated UK10K patients can not be shared because of the UK10K data access policy. 
The annotated UK10K exome files can be shared upon request if the requester can prove to have authorized access to the data. 

### Python dependencies

Python>=3.7
The following scripts have only been tested on unix environments (MacOS and linux)
and we do not offer support for windows users.

We recommend installing the library dependencies inside a conda 
environment by following these steps:

```
conda create --name <env> python=3.7
conda activate <env>

conda install -c conda-forge graph-tool

pip install -r requirements.txt
```


## 2) Running the code

The different types of use of the scripts can be accessed from the `hop.py` script.

We provide for 3 main usages: 
- Running the HOP predictor to predict the cross-validation exomes 
- Running the HOP predictor to prioritize variants
in the independent set exomes
- Running the HOP predictor on any one exome generated from 1KGP and OLIDA combinations
with user-defined seeds.

These three different usages can be perfomed via:

`python hop.py <crossval | indep | predict> `

Different arguments can then be provided depending on the action.

### Common arguments for all actions:

- `-p` `--patient_name` (required): ID of the template patient (it has to be the ID for which the 
annotations are provided)
- `-rf` `--results_folder` (required): The folder in which you want to save the results.
It will be created in the `data` folder. 
- `-r` `--restart` (optional): The value of the restart parameter for the RWR algorithm
  (default = 0.3), value has to be between 0 and 1. 
- `-n` `--nb_top` (optional): Number of combinations to output as top ranked combinations. Default is 100. 
- `-wr` `--write_raw` (optional): Use the flag if you want the raw results to be written, i.e. all score values for all combinatons in the exome
should be written. Defaults to False (see Important Note below)

Important note: With the "write_raw" option, the predictor outputs a file with all generated combinations and
their pathogenicity, disease-relevance and final scores. For some exomes, this can generate extremely large
files even if they are gzipped. If you are running this code on your own laptop or have little space available, 
it might be preferable to not use this option. 

### Reproducing paper results

#### Cross-validation results
The action `crossval` will, for one exome template, divide the training set in 10 folds, train 10 VarCoPP2.0 models, 
and sequentially insert an OLIDA combination in the exome template, predict the exome with the relevant
VarCoPP2.0 model and HOP scores. 

Specific arguments:
- `-st` `--seed_type` (optional): Which type of seeds should be use by HOP. Options
are "All", "HPO", "Panel" and "HPO+panel". Defaults to "All"
- `-o` `--all_operators` (optional): Boolean option for whether you want the different types of operators
  (MAX, MIN, AVGE, MULT) to be computed for the computation of the FinalScore. Defaults to False.
- `-fs` `--final_score` (optional): In the case when "all" is used for the seed type and/or if the "all_operators" 
option is selected, this parameter allows to define the seed and operator used for the ranking of the combinations
in the top N output. Options can be specified as <seed_type>_<operator>. Defaults to HPO+Panel_AVGE.

Example:

`python hop.py crossval -p <patient_id> -rf results_crossval -n 500 -st All`

This will create 301 files in the `data/results_crossval` folder which will contain the top
500 prioritized combinations of each of the synthetic exomes generated from the exome template and the 
OLIDA combinations belonging to the independent test set. 

In addition, the training set will have been splitted in 10 folds saved in the "crossval_folds" folders
and 10 VarCoPP2.0 models, corresponding to each of the folds will have been saved in the "varcopp_data" folder.
This allows to run the same cross-validation on the same exome templates.

#### Independent set results 
The action `indep` will, for one template exome, insert sequentially the OLIDA
 combinations from the independent test set in the exome, rank combinations using HOP 
with the final VarCoPP2.0 model.

- `-st` `--seed_type`: Which type of seeds should be use by HOP. Options
are are "All", "HPO", "Panel" and "HPO+Panel".
- `-fs` `--final_score` (optional): In the case when "all" is used for the seed type and/or if the "all_operators" 
option is selected, this parameter allows to define the seed and operator used for the ranking of the combinations
in the top N output. Options can be specify as <seed_type>_<operator>. Defaults to HPO+Panel_AVGE.

Example:

`python hop.py indep -p <patient_id> -rf results_independent -n 500 -st All`

This will output 119 files in the `data/results_independent` folder which will contain the top
500 prioritized combinations of each of the synthetic exomes generated from the exome template and the 
OLIDA combinations belonging to the independent test set. 

### Prioritizing a single exome

The `predict` action will run HOP on a single exome created by inserting one OLIDA
combination in one template exome, and using the final VarCoPP2.0 predictor to compute the 
Pathogenicity Score.

In addition to the main arguments you need to specify:

- `-c` `--combination` (required): ID of the OLIDA combination to insert in the 
template exome
- `-s` `--seeds` (required): List of seeds to use for the RWR, separated by a comma. Should be a
list of HPO id terms or Ensembl gene IDs. If the seeds are not found in the KG, 
the output ranking will be based on the Pathogenicity Score alone.

Example:

`python hop.py predict -p <patient_id> -rf results_single -n 500 -c OLI002 -s HP:0002092`

This will output one file in the `results/results_single` folder with the top 500 ranked combinations
when the OLI002 is inserted in the exome of the patient with template <patient_id>. 

## 4) Contributors

Main contributor to the code: Barbara Gravel

People who also contributed to the development of the code:
- Alexandre Renaux

## 5) License

This code is licensed under the MIT License (https://choosealicense.com/licenses/mit/).

