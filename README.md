# prepDyn: Preprocessing sequences for dynamic homology

A collection of Python scripts to facilitate the preprocessing of input sequences for dynamic homology. 

In dynamic homology, data should be preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events to avoid grouping from artifacts. However, previous empirical studies using POY/PhyG manually preprocessed data with varying approaches. Here we present **prepDyn**, a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG.

Copyright (C) Daniel Y. M. Nakamura 2025

## Installation

The two dependencies that should be installed beforehand by the user are:
- Python v. 3.10.9 (or newer), including *argparse*, *ast*, *csv*, *importlib*, *re*, *StringIO*, *subprocess*, *sys*, *tempfile*, and *time*, which are usually part of recent versions of Python.
- MAFFT v. 7.5.2 (or newer), installed in $PATH as 'mafft'.

Other dependencies are Python modules that will be automatically installed by prepDyn (if already installed, they will only be loaded):
- Bio v. 1.73 (or newer), including *AlignIO*, *Entrez*, *SeqIO*, *Align*, *Seq*, and *SeqRecord*.
- matplotlib v. 3.7.0 (or newer)
- numpy v. 1.23.5 (or newer)
- termolor

## Introduction

The four steps are (1) data collection from GenBank, (2) trimming, (3) identification of missing data, and (4) successive partitioning.

## Usage
prepDyn is organized in three Python files:
- prepDyn_auxiliary.py: script containing all auxiliary Python functions required by the other scripts.
- GB2MSA.py: script to download sequences from GenBank and identify internal missing data.
- prepDyn.py: script integrating the pipeline.

The following examples are designed for users with little experience on Unix. If you have questions, send a message using **GitHub issues**.

### Example 1: A single alignment

### Example 2: Multiple alignments

### Example 3: Appending new sequences

Given a CSV file called *input.csv*, whose first column is called *Terminals* and the other columns are the names of genes (and each cell contain the correspondent GenBank accession number), the following command will download the sequences, align them with MAFFT and trim orphan nucleotides of length >10 bp. In addition, files containing the names of the terminals (useful for control of taxon sampling in POY/PhyG) and the wall-clock time will be reported. 

```
python GB2MSA.py --input_file input.csv --output_prefix output --orphan_threshold 10 --delimiter , --write_names --log
```

If more than one GenBank accession number is specified in the same cell refering to non-overlapping fragments of the same gene (e.g. MT893619/MT895696), the space between them is identified as missing data (?).

## Citation
