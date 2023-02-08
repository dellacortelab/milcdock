# MILCDock

<img src="https://github.com/dellacortelab/milcdock/blob/main/results/figures/reni_dude-lit-pcba_ensemble_hist.png" alt="drawing" width="500"/>

MILCDock is an consensus docking neural network for virtual screening and drug discovery. This code is intended to be used for structure-based virtual screening to sort databases of thousands to millions of ligands according to their likelihood of activity against a target protein. Using machine learning to combine the predictions of multiple docking tools significantly improves the performance of docking programs on this task (see [paper](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00705)). 

Try out the MILCDock model here:
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dellacortelab/milcdock/blob/main/examples/milc_dock_eval.ipynb)

# Quickstart

The only requirement to run MILCDock is [Docker](https://www.docker.com/), which is used to manage dependencies. Here is one example of the full docking workflow for a small example batch of ligands (NOTE: Make sure to accept the terms of the PLANTS academic license at the [PLANTS website](http://www.tcd.uni-konstanz.de/plants_download/))

Running consensus docking involves three steps:
1. Run docking simulations using all 5 docking programs.
2. Generate machine learning input vectors.
3. Use machine learning to accurately predict the probability each ligand is active.

```
# Create a docker image with the necessary dependencies
docker build -t milcdock -f dependencies/Dockerfile dependencies
# Run a docker container
docker run -dit --gpus all --name docking --rm -v $(pwd):/code milcdock
docker exec -it docking /bin/bash
# (Inside docker container): Download ZINC ligands
cd code
python3 download_zinc.py
# Run docking with all 5 tools on a directory of ligands
./dock.sh -r ./data/0_ligand_receptor/receptor/ace.pdb -c ./data/0_ligand_receptor/receptor/ace_crystal_ligand.mol2 -l ./data/0_ligand_receptor/ligand/zinc_ligands -o ./data/1_docked_output/ --vina --ledock --rdock --plants --autodock
# Generate ML inputs from docked poses
python3 generate_ml_input.py --docking_output_dir ./data/1_docked_output/ace --ml_pkl_path ./data/3_ml_inputs/ml_inputs.pkl
# Run MILCDock neural network to predict inactives and actives (0-1)
python3 milcdock_predict.py --model_load_path ./results/saved_models/mlp/dude-lit-pcba.pt --ml_pkl_path ./data/3_ml_inputs/ml_inputs.pkl --score_path ./data/4_milcdock_scores/scores.pkl
```
Scores are stored in a dataframe in a column named after the model, along with scores from individual docking tools.

# More details on each step
Read this if you want to scale up docking to large databases.

## Set up Docker environment
The intended workflow for dependency management is [Docker](https://www.docker.com/). This will prepackage the correct versions of all necessary software.
```
# Create a docker image with the necessary dependencies (NOTE: Make sure to accept the terms of the PLANTS academic license at the PLANTS website: http://www.tcd.uni-konstanz.de/plants_download/)
docker build -t milcdock -f dependencies/Dockerfile dependencies
# Run a docker container
docker run -dit --gpus all --name docking --rm -v $(pwd):/code milcdock
docker exec -it docking /bin/bash
```

## Prepare ligands and receptor
To perform docking for virtual screening, you need a protein receptor (in PDB format), ligands to dock (in mol2 format), and either a crystal ligand (i.e. a ligand located in the desired binding site, in mol2 format, RECOMMENDED) or the coordinates and size of the binding site. An example receptor, crystal ligand, and ligand download script can be found in `./docking/data/0_ligand_receptor`. Run this command to download the example ligands:
```
python3 download_zinc.py
```
To get additional ligands from the ZINC database, go to the [ZINC website](https://zinc20.docking.org/tranches/home/), select the desired subset, and click the download button. Select mol2 + wget, and click download. This will download a script to download the desired subset of files. Place that script in ./data/0_ligand_receptor/ligand, and run `python3 download_zinc.py`.

## Docking example commands:
After dependencies are set up and receptor and ligands are available, docking can be performed with the dock.sh script.

### Dock with all 5 tools using a crystal ligand. 
NOTE: Docking with all 5 tools at once is allowed, but since the code performs docking sequentially (i.e. one tool at a time), it will take many, many days of computing, especially with more than a couple hundred ligands. It is therefore recommended to run each tool in a separate process on a separate CPU node for faster results.
```
path_to_script/dock.sh 
-r path_to_receptor/receptor_file.pdb 
-l path_to_ligands/all_ligands_directory/ 
-c path_to_crystal_ligand/crystal_ligand.mol2 
-o path_to_output_directory/
--vina --ledock --rdock --plants --autodock
```
### Dock with 2 tools using coordinates and size (10 x 10 x 10 angstrom box centered at (1, 2, 3))
```
path_to_script/dock.sh 
-r path_to_receptor/receptor_file.pdb 
-l path_to_ligands/all_ligands_directory/ 
-o path_to_output_directory/
--center 1 2 3
--size 10 10 10
--ledock --plants
```
### Dock with a single tool using a single ligand (only useful for testing workflow)
```
path_to_script/dock.sh
-r path_to_receptor/receptor_file.pdb
-l path_to_ligands/ligands_directory/ligand.mol2
-c path_to_crystal_ligand/crystal_ligand.mol2 
-o path_to_output_directory/
--rdock
```
Running dock.sh will output docking files in the following directory structure: `path_to_output_directory/receptor_file/docking_tool`

## Extract docking statistics

Docking statistics are extracted from docking files by running generate_ml_input.py.
```
python3 generate_ml_input.py --docking_output_dir path_to_docked_output/receptor/ --ml_pkl_path path_to_ml_inputs/ml_inputs.pkl
```
After running generate_ml_input.py, a fully prepared pkl should be created, ready for use in the ML process.

## Run MILCDock

Make predictions on unseen inputs
```
# Test the base MILCDock model
python3 milcdock_predict.py --model_load_path ./results/saved_models/mlp/dude-lit-pcba.pt --ml_pkl_path ./data/3_ml_inputs/ml_inputs.pkl --score_path ./data/4_milcdock_scores/scores.pkl
# Test the MLICDock ensemble model
python3 milcdock_predict.py --model_load_path ./results/saved_models/ensemble_mlp/dude-lit-pcba --ml_pkl_path ./data/3_ml_inputs/ml_inputs.pkl --score_path ./data/4_milcdock_scores/scores.pkl --model_name_suffix ensemble
```
Make predictions on inputs and evaluate according to 'active'/'inactive' labels
```
python3 milcdock_predict.py --model_load_path ./results/saved_models/mlp/dude-lit-pcba.pt --dataset dude-mini --score_path ./data/4_milcdock_scores/scores.pkl --report_metrics
```
Train a model from scratch
```
python3 milcdock_train.py --model_save_path ./results/models/dude-lit-pcba.pt --run_name dude-lit-pcba --train_dataset dude-lit-pcba --val_dataset dude-lit-pcba --device 0
```
- Fine-tune a model
```
python3 milcdock_train.py --model_load_path ./results/saved_models/mlp/dude-lit-pcba.pt --model_save_path ./results/models/dude-lit-pcba-finetune.pt --run_name dude-lit-pcba-finetune --train_dataset lit-pcba --val_dataset lit-pcba dude --device 0
```

# Known issues:

runrDock.sh
- parialcharge obabel - the obabel commands have a --partialcharge gasteiger argument. This failed on one target (pgh2) for unclear reasons. Also, this step is skipped if there is already a mol2 receptor file. This can happen if the user ran PLANTS docking before rdock, and since runPlants.sh doesn't include this --partialcharge argument, this may introduce some inconsistencies.

# Citation
```
@article{doi:10.1021/acs.jcim.2c00705,
author = {Morris, Connor J. and Stern, Jacob A. and Stark, Brenden and Christopherson, Max and Della Corte, Dennis},
title = {MILCDock: Machine Learning Enhanced Consensus Docking for Virtual Screening in Drug Discovery},
journal = {Journal of Chemical Information and Modeling},
volume = {62},
number = {22},
pages = {5342-5350},
year = {2022},
doi = {10.1021/acs.jcim.2c00705},
    note ={PMID: 36342217},
URL = {https://doi.org/10.1021/acs.jcim.2c00705}
}
```
