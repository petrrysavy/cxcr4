# Hierarchical Bandits for Protein Engineering with SpeedPPI

This repository contains the code needed to reproduce the experiments in the paper.

> **Hierarchical Bandits for Protein Engineering with SpeedPPI**
> Petr Ryšavý, Adéla Drahokoupilová, Mark Osborn, Filip Železný
> To be published.

Here, we combine a hierarchical bandit approach with the SpeedPPI to modify
a protein. We test our approach on the CXCR4 protein, which serves as an
entry point of the HIV virus, and its primary function is represented via
SDF-1 binding. Our task is, therefore, to mutate its sequence so that the
SDF-1 binding is preserved, and the binding to HIV GP120 is weakened.

---

## Installation

1. First, set up the SpeedPPI algorithm in the `SpeedPPI` directory. Follow the guide at (https://github.com/patrickbryant1/SpeedPPI).
 Keep the `predict_single.sh` file provided in the `SpeedPPI` directory and download the UniRef30 (https://colabfold.mmseqs.com/) in fthe older
 `./SpeedPPI/data/uniref30_2023/UniRef30_2023_02`.
2. Install necessary Python packages
 ```bash
 pip3 install numpy matplotlib pandas ast scikit-learn jinja2 scipy biopython
 ```
3. Set the base directory in `Settings.py:66` to be the root of this repository.
4. (optional) Replace the CXCR4, GP120, and SDF-1 sequences with your sequences of interest in the `Settings.py` file.

## Usage

1. Run the following command to search the state space. The predictions are added to the existing ones. Alternatively, you can use the `run-slurm.sh` script.
 ```bash
 python3 main.py
 ``` 
2. Plot the results using
 ```bash
 python3 experiments.py
 ```

