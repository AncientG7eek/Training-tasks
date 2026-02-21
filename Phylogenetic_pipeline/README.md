Description:
This repository contains a pipeline conducting a phylogenetic analysis from list of organism names in txt file to genome trees of these organisms.
Use configs/config.yaml to adjust the parameters of analysis (together with input file). 

Usage:

snakemake is needed to run the scripts. 

> sudo apt install snakemake

environment.yml attached in files allows for recreation of environment.

> conda env create -f environment.yml
> conda activate com_gen

Then to begin the analysis, se the config/config.yaml up to your needs and run:

> snakemake -p

The analysis will continue automatically and save the resulting genome trees in result dir.

