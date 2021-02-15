<!---
Copyright 2021 Yong-beom Kwon and Ju-yong Lee. All rights reserved.
-->

<h1 align="center">
<p>MolFinder
</h1>

<h3 align="center">
<p>an evolutionary algorithm for the global optimization of molecular properties and the extensive exploration of chemical space using SMILES
</h3>

MolFinder finds diverse molecules with desired properties efficiently without any training and a large molecular database. Also, it does not require a lot of computing resources. This repository contains MolFinder's code and results. MolFinder consists of a simple executable file, `molfinder`, and crossover and mutation process file `ModSMI.py`.

## Getting Started
### Prerequisites
* `python` == 3.7.*
* `rdkit` >= 2019.09.3.0+
* `pandas`
* `numpy`
* `matplotlib`

### Installation instructions
We recommand to use anaconda3 virtual-env; simplest way
```
conda create -n molfinder_venv python=3.7
conda activate molfinder_venv
conda install -c rdkit rdkit
conda install pandas numpy
```

## Quickstart
### 0. Prepare dataset

### 1. Run MolFinder algorithm
```
molfinder -i some_data
```

## Parameters of MolFinder
* `-i, --input`: (`str`) SMILES file (csv format) used by the model
* `-r, --random-seed`: (`int`, None) Determines the random number that selects the initial molecules
  
* `--bank-size`: (`int`, 100) Bank size used in the algorithm 
* `--seed-size`: (`int`, 60) The number of parent molecules used to generate child molecules
  
* `-dist, --dist-coef`: (`float`, 0.90) Adjust the $D_{avg}$ value

* `--max-round`: (`int`, 150) The maximum number of round
* `-cvg, --convergent-round`: (`int`, 150) Determines how many rounds the Dcut will converge
  
* `-c, --coefficient`: (`float`, 0.96) Coefficient of objective function
* `--target`: (`SMILES: str`, None) Target molecule 

* `-fp, --fp-method`: (`str`, rdkit) Fingerprint method; Morgan or RDKit (default)

* `-v, --verbosity`: Print RDKit Error message.

* Parameters of paper results
```shell
molfinder -r 12345678 --bank-size 1000 --seed-size 600 -dist 0.90 -c 0.994
```

## Set objective fucntion


## References

* paper: In progress..

