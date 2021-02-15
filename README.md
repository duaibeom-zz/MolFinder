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
Prepare a SMILES file in CSV format. **The first column must be SMILES.**  
(Contains headers, It doesn't have to be just only SMILES.)
```
SMILES
CCN(CC)CCN(C(=O)c1ccc(CCC(F)(F)F)cc1)[C@H]1CCS(=O)(=O)C1
C[C@@H](C(=O)N(C)C)N1[C@H]2CC[C@H]1CC(NC(=O)C1C(C)(C)C1(C)C)C2
CCOc1cc(N2C[C@@H]3C(NC(=O)c4ccn(C)n4)[C@H]3C2)ncn1
O=C(NC[C@@H](CO)Cc1cccnc1)c1ccnc(OC2CCC2)c1
O=C([C@@H]1C[C@H]1c1cccnc1)N1CCC(O)(CNCc2ccccn2)CC1
O=C(CCN1C(=O)[C@H]2CCCC[C@@H]2C1=O)NC1CCN(CC(F)(F)F)CC1
...
```

### 1. Run MolFinder algorithm
```
./molfinder -i sample.csv --max-round 5
```

## Parameters of MolFinder
* `-i, --input`: (`str`) SMILES file (csv format) used by the model
* `-r, --random-seed`: (`int`, None) Determines the random number that selects the initial molecules
  
* `--bank-size`: (`int`, 100) Bank size used in the algorithm 
* `--seed-size`: (`int`, 60) The number of parent molecules used to generate child molecules
  
* `-dist, --dist-coef`: (`float`, 0.90) Adjust the $D_{avg}$ value

* `--max-round`: (`int`, 150) The maximum number of round
* `-cvg, --convergent-round`: (`int`, 150) Determines how many rounds the Dcut will converge
  
* `-c, --coefficient`: (`float`, 0.9) Coefficient of objective function
* `--target`: (`SMILES: str`, None) Target molecule 

* `-fp, --fp-method`: (`str`, rdkit) Fingerprint method; Morgan or RDKit (default)

* `-v, --verbosity`: Print RDKit Error message.

```shell
# Parameters of paper results
molfinder -r 12345678 --bank-size 1000 --seed-size 600 -dist 0.90 -c 0.994 -i sample.csv
```

## Set objective fucntion
1. Find `@@FEATURES` in `molfinder`, Set your features.
2. Find `@@REWARD` in `molfinder`, Modifiy your objective function.

## Component of MolFinder
```
MolFinder Algorithm
├── molfinder
│   ├── cal_avg_dist
│   ├── prepare_seed
│   ├── prepare_child
|   │   └── Crossover and Mutations from ModSMI
│   └── update_bank
└── ModSMI.py
    ├── tight_rm_branch
    ├── prepare_rigid_crossover
    ├── replace_atom
    ├── add_atom
    └── delete_atom
```

## References

* paper: In progress..

---
This is my first code and was mainly written in February 2020. There are many drawbacks, but I keep learning and trying.
