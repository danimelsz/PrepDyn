# prepDyn

A collection of Python scripts to facilitate the preprocessing of input sequences to dynamic homology. 

In dynamic homology, data are typically preprocessed to distinguish differences in sequence length resulting from missing data or insertion-deletion events. However, previous empirical studies using POY/PhyG manually preprocessed missing data with varying approaches. Here we demonstrate that coding missing data with dashes (–) or IUPAC Ns increase tree costs and are biologically inappropriate. Although inserting pound signs (#) around blocks of missing data has been a common solution, it reduces the severity of homology tests and precludes the discovery of the optimal tree. Therefore, missing data should be coded with question marks (?) to minimize tree costs, whereas pound signs should be inserted only on highly conserved regions. To balance time complexity and severity of homology test, we propose a heuristic to successively partition data. All procedures are implemented in a collection of Python scripts to facilitate the preprocessing of input sequences to POY/PhyG. 

## Installation

## Usage

## Citation
