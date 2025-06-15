# prepDyn: Preprocessing sequences for dynamic homology

[![language](https://img.shields.io/badge/language-python-green?style=flat&logo=python&logoColor=green)](https://www.python.org)
[![author](https://img.shields.io/badge/author-DYM_Nakamura-green?logo=googlescholar&logoColor=green)](https://scholar.google.com/citations?user=c0W8Cm8AAAAJ&hl=en)
[![license](https://img.shields.io/badge/license-GPL_v3-green?logo=gnu&logoColor=green)](https://www.gnu.org/licenses/gpl-3.0.html)

A collection of Python scripts to facilitate the preprocessing of input sequences for dynamic homology. 

In dynamic homology, data should be preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events to avoid grouping from artifacts. However, previous empirical studies using POY/PhyG manually preprocessed data with varying approaches. Here we present **prepDyn**, a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG.

Copyright (C) Daniel Y. M. Nakamura 2025

## Installation

The two dependencies that should be installed beforehand by the user are:
- Python v. 3.10.9 (or newer), including *argparse*, *ast*, *csv*, *importlib*, *re*, *StringIO*, *subprocess*, *sys*, *tempfile*, and *time*, which are usually part of recent versions of Python.
- MAFFT v. 7.5.2 (or newer), installed in $PATH as 'mafft'.

```
conda create -n new_env python=3.10 --yes
conda install bioconda::mafft
```

Other dependencies are Python modules that will be automatically installed by **prepDyn** when you run it for the first time:
- Bio v. 1.73 (or newer), including *AlignIO*, *Entrez*, *SeqIO*, *Align*, *Seq*, and *SeqRecord*.
- matplotlib v. 3.7.0 (or newer)
- numpy v. 1.23.5 (or newer)
- termcolor

If the  modules are not installed automatically, try:
```
conda install conda-forge::biopython
conda install conda-forge::matplotlib
conda install anaconda::numpy
conda install conda-forge::termcolor
```

Finally, clone the **prepDyn** repository using the command:
```
git clone https://github.com/danimelsz/PrepDyn.git
```

## Introduction

**prepDyn** comprises four steps: (1) data collection from GenBank, (2) trimming, (3) identification of missing data, and (4) successive partitioning.

## Usage
**prepDyn** is organized in three Python files:
- prepDyn.py: main script integrating the pipeline.
- GB2MSA.py: script to download sequences from GenBank and identify internal missing data.
- prepDyn_auxiliary.py: script containing all auxiliary Python functions required by the other scripts.

The following examples are designed for users with little experience on Unix. If you have questions, send a message using **GitHub issues**.

### Example 1: Basic

The basic use of **prepDyn** is running all four steps using a single command. Given an input CSV, whose first column is called *Terminals* and the other columns are the names of genes (each cell containing the correspondent GenBank accession number), the following command will download sequences, trim invariants and orphan nucleotides <10 bp in terminal positions, and identify missing data as *?* (all differences in sequence length in terminal positions are missing data). The log reports the run time.

```
python prepDyn.py --GB_input input.csv --output_file out --del_inv --orphan_method semi --orphan_threshold 10 --partitioning_round 0 --log
```

In the CSV file, if more than one GenBank accession number is specified in the same cell refering to non-overlapping fragments of the same gene (e.g. MT893619/MT895696), the space between them is automatically identified as internal missing data (?).

We specified *--paritioning_round 0*, which means that partitioning was not performed. As a heuristic, we recommend testing the impact of adding pound signs to the tree optimality scores using a successive partitioning strategy. For instance, if you specify *--partitioning_round 1*, the largest block(s) of contiguous invariants will be partitioned.

```
python prepDyn.py --input_file data.fasta --output_file out1 --partitioning_round 1 --log
```

This process can continue until tree costs reported by POY/PhyG remain stationary (e.g. *--partitioning_round 2* inserts pound signs in the 1- and 2-largest block(s) of contiguous invariants).

```
python prepDyn.py --input_file data.fasta --output_file out2 --partitioning_round 2 --log
```

### Example 2: GB2MSA + prepDyn

Suppose you want to download sequences and preprocess them using different commands. Given a CSV file called *input.csv*, the following command will download the sequences and align them with MAFFT. In addition, files containing the names of the terminals (useful for control of taxon sampling in POY/PhyG) and the run time will be reported. 

```
python GB2MSA.py --input_file input.csv --output_prefix output --delimiter , --write_names --log
```

Now, you can run **prepDyn**:

```
python prepDyn.py --input_file output.fasta --del_inv --orphan_method semi --partitioning_round 0 --log
```

### Example 3: Multiple alignments

Suppose you have a phylogenomic dataset with hundreds of gene alignmens in the directory *./data/*. You can preprocess all gene alignments in FASTA format using a single command:

```
python prepDyn.py --input_file ./data/ --input_format fasta --output_prefix output --del_inv --orphan_method semi --log
```

### Example 4: Appending new sequences

If you have newly generated sequences (unavailable in GenBank) that you want to append to a prealigned FASTA alignment, run:

```
AAAAAAAA
```

### Example 5: Ancient DNA

Suppose you have a dataset with ancient DNA sequences from the sample *Dendropsophus_tritaeniatus_MZUSP73973*. The IUPAC Ns present in sequences are ambiguous positions resulting from low coverage death in DNA read mapping. It is unknown if this positions actually correspond to nucleotides N or to indels (-). As such, the IUPAC Ns of ancient DNA sequences can be replaced with question marks using the following command:

```
AAAAAAAA
```


## Citation
