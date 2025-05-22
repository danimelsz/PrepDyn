#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# auxiliary.py
# prepDyn - Copyright (C) 2025
# Daniel Y. M. Nakamura
# GNU General Public Licence version 3.0
# Contact: dani_ymn@outlook.com

COPYRIGHT=""""
GB2MSA.py is part of prepDyn.

prepDyn - Data preprocessing for dynamic homology
Copyright (C) 2025 - Daniel Y. M. Nakamura
    prepDyn comes with ABSOLUTELY NO WARRANTY!
    This is a free software, and you are welcome to
    redistribute them under certain conditions
    For additional information, please visit:
    http://opensource.org/licenses/GPL-3.0
"""
print(COPYRIGHT)

###########
# MODULES #
###########

# The programs pip and mafft should be installed beforehand.
# The Python modules can be automatically installed using the
# following script. If already installed, it will load them.

import subprocess
import sys
import importlib

def install_and_import(package_name, import_path, from_list=None, alias=None):
    try:
        if from_list:
            module = importlib.import_module(import_path)
            for name in from_list:
                globals()[name] = getattr(module, name)
        else:
            globals()[alias or import_path] = importlib.import_module(import_path)
    except ModuleNotFoundError:
        print(f"Installing: {package_name}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package_name])
        # Try importing again after installation
        if from_list:
            module = importlib.import_module(import_path)
            for name in from_list:
                globals()[name] = getattr(module, name)
        else:
            globals()[alias or import_path] = importlib.import_module(import_path)

# Biopython modules
install_and_import("biopython", "Bio.AlignIO", from_list=["read", "write"])
install_and_import("biopython", "Bio.Entrez", from_list=["efetch", "email"])
install_and_import("biopython", "Bio.SeqIO", from_list=["parse", "write"])
install_and_import("biopython", "Bio.Align", from_list=["MultipleSeqAlignment"])
install_and_import("biopython", "Bio.Seq", from_list=["Seq"])
install_and_import("biopython", "Bio.SeqRecord", from_list=["SeqRecord"])

# Other useful packages
install_and_import("matplotlib", "matplotlib.pyplot", alias="plt")
install_and_import("numpy", "numpy", alias="np")
install_and_import("termcolor", "termcolor", from_list=["colored"])

# Standard libraries (should not need installation)
import argparse
import ast
import csv
from io import StringIO
import os
import pathlib
import re
import subprocess
import tempfile
import time

# prepDyn libraries
from prepDyn_auxiliary import GB2MSA
from Bio import SeqIO
from Bio.Seq import Seq

########
# MAIN #
########

def main():
    parser = argparse.ArgumentParser(
        description="GB2MSA downloads sequences from GenBank and performs multiple sequence alignment using MAFFT. If >1 non-overlapping fragments of the same gene is specified (e.g. MT893619/MT895696), the space between them is identified as missing data (?)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""\
Example:  python GB2MSA.py --input_file input.csv --output_prefix myoutput --delimiter , --write_names --log --orphan_threshold 10
""")


    parser.add_argument("--input_file", type=str, required=True, help="Path to CSV/TSV file with GenBank accession numbers.")
    parser.add_argument("--output_prefix", type=str, required=True, help="Path (including prefix, if desirable) for output FASTA files.")
    parser.add_argument("--delimiter", default=",", type=str, required=False, help="Delimiter used in the input file (default: ',').")
    parser.add_argument("--write_names", action="store_true",  required=False, help="If set, write sequence names in a separate file, which can be used as input data in POY/PhyG to select taxon sample (default: False).")
    parser.add_argument("--log", action="store_true",  required=False, help="If set, write wall and CPU time to a log file (default: False).")
    parser.add_argument("--orphan_threshold", type=int, default=10,  required=False, help="Threshold to clean orphan nucleotides (default: 10).")

    args = parser.parse_args()

    GB2MSA(
        input_file=args.input_file,
        output_prefix=args.output_prefix,
        delimiter=args.delimiter,
        write_names=args.write_names,
        log=args.log,
        orphan_threshold=args.orphan_threshold)

# Only run main() directly from command line, not imported as a module
if __name__ == "__main__":
    main()