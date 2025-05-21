# prepDyn: Preprocessing sequences for dynamic homology

A collection of Python scripts to facilitate the preprocessing of input sequences to dynamic homology. 

In dynamic homology, data should be preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events to avoid grouping from artifacts. However, previous empirical studies using POY/PhyG manually preprocessed data with varying approaches. Here we present **prepDyn**, a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG. The four steps are (1) data collection from GenBank, (2) trimming, (3) identification of missing data, and (4) successive partitioning.

Copyright (C) Daniel Y. M. Nakamura 2025

## Installation

The two dependencies that should be installed beforehand by the user are:
- Python v. 3.10.9 (including *argparse*, *ast*, *csv*, *importlib*, *re*, *StringIO*, *subprocess*, *sys*, *tempfile*, and *time*, which are usually part of recent versions of Python).
- MAFFT v. 7.5.2 (or later), installed in $PATH as 'mafft'.

Other dependencies are the following Python modules that will be automatically installed by prepDyn (if already installed, they will only be loaded):
- Bio (including *AlignIO*, *Entrez*, *SeqIO*, *Align*)
- matplotlib
- numpy
- termolor

## Usage

If you have questions, send a message using **GitHub issues**.

## Citation
