#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# auxiliary.py
# prepDyn - Copyright (C) 2025
# Daniel Y. M. Nakamura
# GNU General Public Licence version 3.0
# Contact: dani_ymn@outlook.com

COPYRIGHT=""""
auxialiary.py is part of prepDyn.

prepDyn - Data preprocessing for dynamic homology
Copyright (C) 2025 - Daniel Y. M. Nakamura
    prepDyn comes with ABSOLUTELY NO WARRANTY!
    This is a free software, and you are welcome to
    redistribute them under certain conditions
    For additional information, please visit:
    http://opensource.org/licenses/GPL-3.0
"""

###########
# MODULES #
###########

# The programs pip and Mafft should be installed beforehand.
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

# Load packages
from Bio import AlignIO, Entrez, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import os
import pathlib
import re
import subprocess
import tempfile
from termcolor import colored
import time

#######################
# AUXILIARY FUNCTIONS #
#######################

# Function 1. Visualize colored alignments
# Define colors for each nucleotide (case insensitive)
def color_nucleotide(nucleotide):
    color_map = {
        "A": "red",
        "T": "blue",
        "G": "green",
        "C": "yellow",
        "-": "white",
        "#": "black",
        "?": "magenta"
    }
    return colored(nucleotide, color_map.get(nucleotide.upper(), "white"))
# Print the alignment
def print_colored_alignment(alignment):
    """
    Prints a DNA alignment with each nucleotide in a different color.

    Parameters:
        alignment (dict): A dictionary with sequence names as keys and DNA sequences as values.
    """
    for name, sequence in alignment.items():
        colored_seq = ''.join(color_nucleotide(n) for n in sequence)
        print(f"{name}: {colored_seq}")

# Function 2. Convert dict to MultipleSeqAlignment
def dict_to_multiple_seq_alignment(seq_dict):
    """
    Converts a dictionary of sequences into a MultipleSeqAlignment object.
    
    Args:
        seq_dict (dict): A dictionary where keys are sequence identifiers and values are sequences (strings).
        
    Returns:
        MultipleSeqAlignment: A Biopython MultipleSeqAlignment object.
    """
    # Create a list of SeqRecord objects from the dictionary
    seq_records = []
    
    for seq_id, seq_str in seq_dict.items():
        # Create a SeqRecord for each sequence
        seq = Seq(seq_str)
        seq_record = SeqRecord(seq, id=seq_id, description="")
        seq_records.append(seq_record)
    
    # Create and return the MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(seq_records)
    return alignment

# Function 3. List lengths of blocks of contiguous gaps in internal and terminal positions
def list_gap_blocks_by_type(alignment, plot_distribution=False):
    """
    Identify all blocks of contiguous gaps in the DNA alignment.
    Classify gap blocks into terminal and internal blocks.
    Optionally, plot the distributions of gap block lengths for terminal and internal blocks.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        plot_distribution (bool): If True, plot the distributions of terminal and internal gap block lengths.
        
    Returns:
        tuple: Two lists:
            - terminal_blocks: A list of gap block lengths at the start or end of sequences.
            - internal_blocks: A list of gap block lengths in the middle of sequences.
    """
    terminal_blocks = []
    internal_blocks = []
    
    # Iterate over each sequence in the alignment
    for seq_id, sequence in alignment.items():
        sequence_length = len(sequence)
        gap_count = 0
        in_gap = False

        # Identify gap blocks and separate terminal vs internal
        for i, nucleotide in enumerate(sequence):
            if nucleotide == '-':  # We are in a gap
                if not in_gap:
                    in_gap = True
                    gap_count = 1  # Start a new gap block
                else:
                    gap_count += 1  # Continue counting the current gap block
            else:  # We encountered a non-gap nucleotide
                if in_gap:
                    # End of a gap block
                    # Check if this gap block is terminal (at start or end of the sequence)
                    if i == sequence_length or (i - gap_count == 0):
                        terminal_blocks.append(gap_count)
                    else:
                        internal_blocks.append(gap_count)
                    in_gap = False
                    gap_count = 0  # Reset gap count for the next block

        # If the sequence ends with a gap block, add it to the appropriate list
        if in_gap:
            if sequence[-1] == '-':  # Last character is a gap
                terminal_blocks.append(gap_count)
            else:
                internal_blocks.append(gap_count)

    # Optionally plot the distributions of terminal and internal blocks
    if plot_distribution:
        plt.figure(figsize=(8, 6))
        # Plotting terminal blocks
        plt.hist(terminal_blocks, bins=20, alpha=0.5, label='Terminal Blocks', color='blue')
        # Plotting internal blocks
        plt.hist(internal_blocks, bins=20, alpha=0.5, label='Internal Blocks', color='red')
        
        plt.xlabel('Gap Block Length')
        plt.ylabel('Frequency')
        plt.title('Distribution of Gap Block Lengths')
        plt.legend(loc='upper right')
        plt.show()

    return terminal_blocks, internal_blocks

# Function 4. Remove underscores in the beggining of file names
def remove_leading_underscores(file_path):
    """
    Remove leading contiguous underscores from the beginning of the file name.
    If a directory is provided, it renames all files in that directory.
    
    Args:
        file_path (str): The full path of the file or directory.
    
    Returns:
        str or None: The new file path with leading underscores removed if a file, 
                     or None if a directory (modifies in place).
    """
    if os.path.isdir(file_path):
        # If it's a directory, rename all files inside it
        for root, dirs, files in os.walk(file_path):
            for file in files:
                old_file_path = os.path.join(root, file)
                new_file_name = file.lstrip('_')
                new_file_path = os.path.join(root, new_file_name)
                
                if old_file_path != new_file_path:
                    os.rename(old_file_path, new_file_path)
        return None  # No return for directories, as the renaming is in place
    else:
        # If it's a single file, rename it
        dir_name, file_name = os.path.split(file_path)
        new_file_name = file_name.lstrip('_')
        new_file_path = os.path.join(dir_name, new_file_name)
        if file_path != new_file_path:
            os.rename(file_path, new_file_path)
        return new_file_path
    
############################################   
# MAIN FUNCTIONS: STEP 1. DATA COLLECTION  #
############################################

def GB2MSA_1(input_file, output_prefix, delimiter=',', write_names=True):
    """
    Downloads GenBank sequences based on accession numbers in a CSV/TSV file and aligns them by gene using 
    MAFFT. If two fragments of the same locus are concatenated with no overlap between them, the space
    between them will be treaed as missing data (15 Ws will flank these blocks of missing data, which will
    be used to track these regions and be replaced with question marks in GB2MSA_2).

    Parameters:
    -----------
    input_file : str
        Path to the CSV or TSV input file. The first column should contain sequence names (sample identifiers).
        The first row should contain gene names starting from the second column. Cells contain GenBank 
        accession numbers (one or more separated by slashes). "NA", empty cells, or dashes are ignored.

    output_prefix : str
        Prefix used for naming intermediate FASTA files and final aligned output files.

    delimiter : str, optional (default=',')
        Delimiter used in the input file (e.g., ',' for CSV or '\t' for TSV).

    write_names : bool, optional (default=True)
        If True, writes a TXT file listing all sequence names (from the first column).

    Returns:
    --------
    aligned_files : list of str
        List of file paths to the MAFFT-aligned FASTA files for each gene.
    """
    # Open the input CSV/TSV file
    with open(input_file, newline='') as file:
        reader = csv.reader(file, delimiter=delimiter)
        rows = list(reader)

    # Replace spaces in sequence names with underscores
    sequence_names = [row[0].replace(" ", "_") for row in rows[1:]]
    
    # Extract gene names from the header row (excluding first column)
    gene_names = rows[0][1:]
    
    # Extract gene accession data for each sequence (excluding first column)
    gene_columns = [row[1:] for row in rows[1:]]

    # If requested, write the sequence names to a text file
    if write_names:
        names_file = f"{output_prefix}_sequence_names.txt"
        with open(names_file, 'w') as nf:
            for name in sequence_names:
                nf.write(f"{name}\n")

    aligned_files = []  # List to store paths of aligned output files

    # Iterate over each gene (i.e., each column after the first)
    for gene_idx, gene_name in enumerate(gene_names):
        fasta_file = f"{output_prefix}_{gene_name}.fasta"  # Name of temporary FASTA file
        
        with open(fasta_file, 'w') as fasta_out:
            # Iterate through each row (sample/sequence)
            for i, seq_name in enumerate(sequence_names):
                cell = gene_columns[i][gene_idx].strip()
                
                # Skip cells with missing data ("NA", empty, or dash)
                if cell.upper() == "NA" or not cell or cell == "-":
                    continue

                # Split accession numbers by '/' and filter out invalid entries
                accessions = [acc for acc in cell.split('/') if acc.upper() != "NA" and acc != "" and acc != "-"]
                sequences = []  # To hold the sequences retrieved from GenBank

                # Fetch each sequence from GenBank
                for acc in accessions:
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
                        seq_record = SeqIO.read(handle, "fasta")
                        handle.close()
                        sequences.append(str(seq_record.seq))  # Store the sequence string
                    except Exception as e:
                        print(f"Error fetching {acc}: {e}")

                # Combine multiple sequences with 'W' delimiters to mark junctions (to be handled later)
                combined_seq = "WWWWWWWWWWWWWWW".join(sequences)
                
                # Write the sequence to the FASTA file with its name as header
                fasta_out.write(f">{seq_name}\n{combined_seq}\n")

        # Define the name for the output alignment file
        aligned_file = f"{output_prefix}_{gene_name}_aligned.fasta"

        # Run MAFFT on the generated FASTA file and save the alignment
        with open(aligned_file, 'w') as aligned_out:
            subprocess.run(["mafft", "--auto", fasta_file], stdout=aligned_out)

        # Append the aligned file path to the result list
        aligned_files.append(aligned_file)

    return aligned_files  # Return list of aligned output file paths

def GB2MSA_2(alignment_file):
    """
    If internal missing data were identified by GB2MSA_1, 15 Ws flank the blocks of missing data. 
    For each sequence in the alignment, GB2MSA_2:
    - Replaces internal blocks of 15 'w' or spaced 'w' (e.g., w-w-w) with question marks.
    - Removes columns with only '?' or '-' in all rows.
    - Replaces dash blocks flanked by a nucleotide and a question mark (e.g., A??--AAC) with 
    question marks. These dashes can be artifacts from spaced 'w' and are missing data.

    Parameters:
    -----------
    alignment_file : str
        Path to the MAFFT-aligned FASTA file to be processed.
    
    Returns:
    --------
    str
        Path to the cleaned alignment file.
    """
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    updated_records = []

    # Step 1: Replace exact block of 15 'w' or 'W' with 15 '?'
    for record in alignment:
        seq = str(record.seq)
        seq_cleaned = seq.replace("w" * 15, "?" * 15).replace("W" * 15, "?" * 15)
        record.seq = Seq(seq_cleaned)
        updated_records.append(record)

    # Step 2: Replace non-contiguous 'w' blocks (e.g., w-w-w) with '?'
    for record in updated_records:
        chars = list(str(record.seq))
        i = 0
        while i < len(chars):
            if chars[i].lower() == 'w':
                count = 1
                indices = [i]
                j = i + 1
                while j < len(chars) and count < 15:
                    if chars[j] == '-':
                        indices.append(j)
                    elif chars[j].lower() == 'w':
                        indices.append(j)
                        count += 1
                    else:
                        break
                    j += 1
                if count == 15:
                    for idx in indices:
                        chars[idx] = '?'
                i = j
            else:
                i += 1
        record.seq = Seq(''.join(chars))

    # Step 3: Remove columns with only '?' or '-' in all rows
    sequences = [list(str(record.seq)) for record in updated_records]
    if len(set(len(seq) for seq in sequences)) > 1:
        raise ValueError("Sequences are not of the same length!")

    valid_columns = []
    for i in range(len(sequences[0])):
        column = [seq[i] for seq in sequences]
        if any(base not in ['?', '-'] for base in column):
            valid_columns.append(i)

    # Rebuild records with cleaned sequences
    cleaned_records = []
    for record in updated_records:
        cleaned_seq = ''.join(str(record.seq)[i] for i in valid_columns)
        record.seq = Seq(cleaned_seq)
        cleaned_records.append(record)

    # Step 4: Replace dash blocks flanked by nucleotide and '?' with '?'
    for record in cleaned_records:
        chars = list(str(record.seq))
        i = 0
        while i < len(chars):
            if chars[i] == '-':
                start = i
                while i < len(chars) and chars[i] == '-':
                    i += 1
                end = i - 1

                # Check flanking characters safely
                left = chars[start - 1] if start > 0 else ''
                right = chars[end + 1] if end + 1 < len(chars) else ''

                if ((left in 'ACGTacgt' and right in '?Nn') or 
                    (left in '?Nn' and right in 'ACGTacgt')):
                    for j in range(start, end + 1):
                        chars[j] = '?'
            else:
                i += 1
        record.seq = Seq(''.join(chars))

    # Step 5: Write to cleaned output
    cleaned_file = alignment_file.replace(".fasta", "_GB2MSA.fasta")
    with open(cleaned_file, "w") as out_handle:
        SeqIO.write(cleaned_records, out_handle, "fasta")

    return cleaned_file

def GB2MSA_3(alignment_dict, orphan_threshold=6, log=False):
    """
    Replaces specific blocks in a DNA alignment with '?':
    - Contiguous gap blocks (length >= orphan_threshold) that are adjacent
      to contiguous nucleotide blocks (length < orphan_threshold), where the
      other side of that nucleotide block touches a '?'

    Parameters:
    alignment_dict (dict): {sequence_name: aligned_sequence}
    orphan_threshold (int): Minimum size to define a valid gap block
    log (bool): If True, print the start and end positions of replaced blocks

    Returns:
    dict: Cleaned alignment with selective replacements
    """
    cleaned_alignment = {}

    for seq_name, original_seq in alignment_dict.items():
        sequence = list(original_seq)
        seq_len = len(sequence)
        updated_seq = ''.join(sequence)

        gap_matches = list(re.finditer(r'-+', updated_seq))
        valid_gap_blocks = [(m.start(), m.end() - 1) for m in gap_matches
                            if (m.end() - m.start()) >= orphan_threshold]

        to_replace = set()

        if log:
            print(f"Sequence: {seq_name}")

        for gap_start, gap_end in valid_gap_blocks:
            replaced = False

            # Check left nucleotide orphan
            left_end = gap_start - 1
            i = left_end
            while i >= 0 and updated_seq[i] not in '-?':
                i -= 1
            left_start = i + 1
            left_len = left_end - left_start + 1

            if left_len > 0 and left_len < orphan_threshold and i >= 0 and updated_seq[i] == '?':
                to_replace.update(range(left_start, left_end + 1))
                to_replace.update(range(gap_start, gap_end + 1))
                replaced = True
                if log:
                    print(f"  Gap block: {gap_start}-{gap_end}")
                    print(f"  Left orphan nucleotide block: {left_start}-{left_end}")

            if not replaced:
                # Check right nucleotide orphan
                right_start = gap_end + 1
                i = right_start
                while i < seq_len and updated_seq[i] not in '-?':
                    i += 1
                right_end = i - 1
                right_len = right_end - right_start + 1

                if right_len > 0 and right_len < orphan_threshold and i < seq_len and updated_seq[i] == '?':
                    to_replace.update(range(right_start, right_end + 1))
                    to_replace.update(range(gap_start, gap_end + 1))
                    if log:
                        print(f"  Gap block: {gap_start}-{gap_end}")
                        print(f"  Right orphan nucleotide block: {right_start}-{right_end}")

        # Replace in sequence
        for i in to_replace:
            sequence[i] = '?'

        cleaned_alignment[seq_name] = ''.join(sequence)

    return cleaned_alignment

def GB2MSA_4(alignment_dict):
    cleaned_alignment = {}

    for seq_name, seq in alignment_dict.items():
        sequence = list(seq)

        # Use regex to find blocks with at least 15 characters of w/W/- combined
        # but must contain at least 15 w or W
        pattern = re.finditer(r'([wW\-]{15,})', ''.join(sequence))

        for match in pattern:
            block = match.group()
            start, end = match.start(), match.end()

            # Count number of w/W in the block
            w_count = sum(1 for c in block if c in 'wW')
            if w_count >= 15:
                for i in range(start, end):
                    sequence[i] = '?'

        cleaned_alignment[seq_name] = ''.join(sequence)

    return cleaned_alignment

####################################
# MAIN FUNCTIONS: STEP 2. TRIMMING #
####################################

def remove_all_gap_columns(alignment):
    """
    Remove columns from the alignment where all terminal sequences have a gap ('-').
    
    Parameters:
        alignment (dict): A dictionary where keys are sequence names and values are sequence strings.
        
    Returns:
        dict: Updated alignment with columns removed where all terminal sequences have a gap.
    """
    # Convert the alignment to a list of sequences
    sequences = list(alignment.values())
    num_sequences = len(sequences)
    seq_length = len(sequences[0])  # Assuming all sequences are the same length

    # Identify columns that need to be removed (where all terminal sequences have a gap)
    columns_to_remove = []
    for col in range(seq_length):
        # Check if all terminal sequences (sp1, sp2, ...) have a gap ('-') in this column
        if all(seq[col] == '-' for seq in sequences):
            columns_to_remove.append(col)
    
    # Remove the identified columns from each sequence
    for seq_name, seq in alignment.items():
        # Create a new sequence with the columns removed
        new_seq = ''.join(seq[col] for col in range(seq_length) if col not in columns_to_remove)
        alignment[seq_name] = new_seq
    
    return alignment

def calculate_orphan_threshold_from_percentile(alignment, percentile=25, log=False, terminal_only=True):
    """
    Calculate orphan threshold based on a specific percentile of gap lengths in the alignment.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        percentile (float): The percentile to use for setting the orphan threshold (e.g., 75 for the 75th percentile).
        log (bool): If True, print the list of gap lengths and the computed orphan threshold.
        terminal_only (bool): If True, only consider gap blocks at the terminal positions (start and end of sequences).
        
    Returns:
        int: Calculated orphan threshold based on the specified percentile.
    """
    gap_lengths = []

    # Loop through each sequence to calculate gap lengths
    for sequence in alignment.values():
        sequence_length = len(sequence)
        gap_count = 0
        in_gap = False

        # If terminal_only is True, we'll check only the first and last blocks of gaps
        if terminal_only:
            # Check for gaps starting at the beginning of the sequence
            if sequence[0] == '-':
                gap_count = 1
                in_gap = True
                for nucleotide in sequence[1:]:
                    if nucleotide == '-':
                        gap_count += 1
                    else:
                        break  # End of the first terminal gap block
                gap_lengths.append(gap_count)  # Add the terminal gap block at the start
                
            # Check for gaps starting at the end of the sequence
            gap_count = 0
            in_gap = False
            if sequence[-1] == '-':
                gap_count = 1
                in_gap = True
                for nucleotide in reversed(sequence[:-1]):
                    if nucleotide == '-':
                        gap_count += 1
                    else:
                        break  # End of the last terminal gap block
                gap_lengths.append(gap_count)  # Add the terminal gap block at the end
        else:
            # Loop through the sequence to find contiguous blocks of gaps (not restricted to terminals)
            for nucleotide in sequence:
                if nucleotide == '-':
                    if not in_gap:
                        in_gap = True  # Start of a new gap block
                        gap_count = 1  # Start counting the length of the new gap block
                    else:
                        gap_count += 1  # Continue counting the gap block length
                else:
                    if in_gap:
                        gap_lengths.append(gap_count)  # End of a gap block, save its length
                        in_gap = False
                        gap_count = 0  # Reset the gap count for the next block

            # If the sequence ends with a gap block, make sure to add the last block's length
            if in_gap:
                gap_lengths.append(gap_count)
        
    # Calculate the specified percentile of gap lengths
    orphan_threshold = int(np.percentile(gap_lengths, percentile))  # Use the given percentile (e.g., 75th percentile)
    
    # Print the list of gap lengths
    if log:
        print("List of gap lengths:", gap_lengths)
        print("Orphan threshold:", orphan_threshold)
        
    return orphan_threshold
    
def delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=False):
    """
    Iteratively eliminates orphan nucleotide blocks from the start and end of each sequence.
    An orphan block is a short contiguous run of nucleotides near the terminal ends separated by many gaps.

    Args:
        alignment (dict): {sequence_id: sequence_string}
        orphan_threshold (int): Max block length and max gap tolerance.
        log_changes (bool): Whether to return a log of changes made.

    Returns:
        dict: Cleaned alignment.
        str (optional): Log of changes made.
    """
    change_log = []

    def find_blocks(seq):
        """Find contiguous non-gap blocks as (start, end) tuples."""
        blocks = []
        i = 0
        while i < len(seq):
            if seq[i] != '-':
                start = i
                while i < len(seq) and seq[i] != '-':
                    i += 1
                end = i
                blocks.append((start, end))
            else:
                i += 1
        return blocks

    alignment_changed = True
    while alignment_changed:
        alignment_changed = False

        for seq_id, sequence in alignment.items():
            seq_list = list(sequence)
            changed = False

            # Iteratively check from left side
            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                first_start, first_end = blocks[0]
                next_start = blocks[1][0]
                gap_count = seq_list[first_end:next_start].count('-')

                if (first_end - first_start < orphan_threshold) and (gap_count > orphan_threshold):
                    deleted = ''.join(seq_list[first_start:first_end])
                    seq_list[first_start:first_end] = ['-'] * (first_end - first_start)
                    changed = True
                    if log_changes:
                        change_log.append(
                            f"{seq_id}: Left block {first_start}-{first_end} deleted ('{deleted}')"
                        )
                else:
                    break

            # Iteratively check from right side
            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                last_start, last_end = blocks[-1]
                prev_end = blocks[-2][1]
                gap_count = seq_list[prev_end:last_start].count('-')

                if (last_end - last_start < orphan_threshold) and (gap_count > orphan_threshold):
                    deleted = ''.join(seq_list[last_start:last_end])
                    seq_list[last_start:last_end] = ['-'] * (last_end - last_start)
                    changed = True
                    if log_changes:
                        change_log.append(
                            f"{seq_id}: Right block {last_start}-{last_end} deleted ('{deleted}')"
                        )
                else:
                    break

            if changed:
                alignment_changed = True
                alignment[seq_id] = ''.join(seq_list)

    if log_changes:
        return alignment, "\n".join(change_log)
    else:
        return alignment

def is_parsimony_non_informative(column):
    """
    Check if a column is non-informative.
    A column is non-informative if:
    - All characters are the same.
    - The characters include '?' and only one other character (e.g., 'A' and '?').

    Parameters:
        column (list): A list of characters in a column (e.g., ['A', 'A', 'A', '?']).
    
    Returns:
        bool: True if the column is non-informative, False if informative.
    """
    unique_characters = set(column)
        
    # Check if the column has only one unique character (parsimony non-informative)
    if len(unique_characters) == 1:
        return True
    
    # Check if the column contains '?' and exactly one other unique character
    if '?' in unique_characters and len(unique_characters) == 2:
        return True
    
    # Otherwise, the column is informative (multiple unique characters without '?')
    return False

def remove_non_informative_positions(alignment, removed_indices=None):
    """
    Remove parsimony non-informative positions from both the start and the end of the alignment.
    
    Parameters:
        alignment (dict): Dictionary where keys are sequence names and values are sequences (strings).
        removed_indices (list, optional): If provided, the function will append the indices of removed columns.
        
    Returns:
        dict: Updated alignment with non-informative positions removed.
    """
    # Convert the alignment to a list of sequences
    sequences = list(alignment.values())
    seq_length = len(sequences[0])

    # Remove non-informative positions from the start
    first_position = 0
    while first_position < seq_length and is_parsimony_non_informative([seq[first_position] for seq in sequences]):
        first_position += 1

    # Remove non-informative positions from the end
    last_position = seq_length - 1
    while last_position >= first_position and is_parsimony_non_informative([seq[last_position] for seq in sequences]):
        last_position -= 1

    # If logging, record the indices of removed columns
    if removed_indices is not None:
        removed_start = list(range(0, first_position))
        removed_end = list(range(last_position + 1, seq_length))
        removed_indices.extend(removed_start + removed_end)

    # Build the updated alignment
    for seq_name, seq in alignment.items():
        alignment[seq_name] = seq[first_position:last_position + 1]

    return alignment

##########################################################
# MAIN FUNCTIONS: STEP 3. IDENTIFICATION OF MISSING DATA #
##########################################################

def replace_terminal_gaps_dict(alignment):
    """
    Replaces terminal gaps ('-') with '?' in a sequence alignment stored as a dictionary,
    while keeping internal gaps as '-'.

    Parameters:
        alignment (dict): Dictionary with sequence names as keys and sequences as values.

    Returns:
        dict: Modified alignment with terminal gaps replaced by '?'.
    """
    # Rename 'modified_alignment' to 'alignment'
    for name, sequence in alignment.items():
        # Find the first and last non-gap characters
        first_non_gap = next((i for i, char in enumerate(sequence) if char != "-"), None)
        last_non_gap = next((i for i, char in enumerate(reversed(sequence), 1) if char != "-"), None)
        
        if first_non_gap is not None and last_non_gap is not None:
            last_non_gap = len(sequence) - last_non_gap  # Adjust reversed index
            
            # Replace terminal gaps with '?' and keep internal gaps as '-'
            new_sequence = (
                "?" * first_non_gap +
                sequence[first_non_gap:last_non_gap + 1] +
                "?" * (len(sequence) - last_non_gap - 1)
            )
            alignment[name] = new_sequence
        else:
            # Handle sequences with only gaps
            alignment[name] = "?" * len(sequence)
    
    return alignment

def replace_dashes_with_question_marks(alignment, internal_column_ranges=None, internal_leaves="all", internal_method="manual", internal_threshold=None):
    """
    Replace dashes with question marks in the specified column ranges for each sequence in the alignment.
    
    Args:
        alignment (dict): A dictionary with sequence IDs as keys and sequences as values.
        internal_column_ranges (list of tuples, optional): A list of tuples where each tuple defines a range of columns 
                                                          (inclusive) to check for dashes. E.g., [(5, 10), (50, 60)].
        internal_leaves (str or list, optional): If "all", replace dashes in all sequences. If a list of sequence IDs 
                                                   is provided, replace dashes in those sequences only.
        internal_method (str, optional): Defines how to replace dashes.
            - "manual": Specify column ranges and terminal sequences to replace dashes.
            - "semi": Replace internal blocks of contiguous dashes larger than the threshold with question marks.
        internal_threshold (int, optional): The threshold for "semi" method. Only internal blocks of contiguous gaps 
                                             larger than this threshold are replaced with question marks.
        
    Returns:
        dict: Updated alignment with dashes replaced by question marks in the specified columns.
    """
    
    # If internal_leaves is a list, only consider those sequences
    if internal_leaves != "all":
        sequences_to_process = set(internal_leaves)
    else:
        sequences_to_process = set(alignment.keys())
    
    # Convert the alignment into a list of sequences for easier indexing
    alignment = {seq_id: list(seq) for seq_id, seq in alignment.items()}  # Convert sequences to lists for mutability

    if internal_method == "manual":
        # Replace dashes with question marks in the specified column ranges
        for seq_id, seq in alignment.items():
            if seq_id not in sequences_to_process:
                continue  # Skip sequences not in the internal_leavesinternal_terminals list
            
            for start, end in internal_column_ranges:
                # Ensure the range is within the bounds of the sequence length
                start = max(0, start)
                end = min(len(seq), end)

                # Replace dashes with question marks within the specified range
                for i in range(start, end + 1):  # +1 because the end is inclusive
                    if seq[i] == '-':
                        seq[i] = '?'
        
    elif internal_method == "semi" and internal_threshold is not None:
        # Replace internal blocks of contiguous dashes larger than the internal_threshold with question marks
        for seq_id, seq in alignment.items():
            if seq_id not in sequences_to_process:
                continue  # Skip sequences not in the internal_leavesinternal_terminals list

            # Identify contiguous blocks of gaps (internal and terminal)
            gap_block_start = None
            for i in range(len(seq)):
                if seq[i] == '-':
                    if gap_block_start is None:
                        gap_block_start = i  # Start of a new gap block
                else:
                    if gap_block_start is not None:
                        # End of a gap block
                        gap_length = i - gap_block_start
                        if gap_length > internal_threshold and gap_block_start != 0 and gap_block_start != len(seq) - gap_length:
                            # It's an internal block larger than threshold, replace with '?'
                            for j in range(gap_block_start, i):
                                seq[j] = '?'
                        gap_block_start = None  # Reset for the next block
            # Check for a gap block at the end of the sequence
            if gap_block_start is not None:
                gap_length = len(seq) - gap_block_start
                if gap_length > internal_threshold and gap_block_start != 0:
                    for j in range(gap_block_start, len(seq)):
                        seq[j] = '?'
        
    # Convert the list back to a string
    alignment = {seq_id: ''.join(seq) for seq_id, seq in alignment.items()}

    return alignment

###################################################
# MAIN FUNCTIONS: STEP 4. SUCCESSIVE PARTITIONING #
###################################################

def add_breaks_terminal(alignment):
    """
    Add # in all instances of terminal gap opening/closure (indicated by ?).
       
    Parameters:
        alignment (dict): Dictionary where keys are sequence names and values are sequences.
    
    Returns:
        dict: Updated alignment with '#' before terminal gap opening and after terminal gap closure.
    """
    # Determine the length of the sequences
    seq_length = len(next(iter(alignment.values())))

    # Initialize a list to keep track of positions that need '#' in all sequences
    hash_positions = [False] * seq_length

    # Iterate through each sequence to find gap regions
    for seq in alignment.values():
        i = 0
        while i < seq_length:
            if seq[i] == '?':
                # Found the start of a gap region
                start = i
                while i < seq_length and seq[i] == '?':
                    i += 1
                end = i
                # Mark the positions for this gap region (avoid marking the first column)
                if start > 0:
                    hash_positions[start] = True
                if end < seq_length:
                    hash_positions[end] = True
            else:
                i += 1

    # Update each sequence with '#' at the identified positions
    for key in alignment:
        new_seq = []
        for i in range(seq_length):
            if hash_positions[i]:
                new_seq.append('#')
            new_seq.append(alignment[key][i])
        # Handle the case where the last position is a gap
        if hash_positions[-1]:
            new_seq.append('#')
        alignment[key] = ''.join(new_seq)

def classify_and_insert_hashtags(alignment, 
                                 partitioning_round=1, 
                                 log_csv_output=False, 
                                 csv_file_path="contiguous_invariant_blocks.csv"):
    # Step 1: Classify columns as invariant or variant
    num_sequences = len(alignment)
    num_columns = len(next(iter(alignment.values())))  # Get the number of columns from one sequence

    column_types = []  # To store the type of each column (invariant/variant)
    contiguous_invariant_blocks = []  # To store lengths and positions of invariant blocks

    for col_idx in range(num_columns):
        column = [seq[col_idx] for seq in alignment.values()]
        unique_values = set(column)
        if len(unique_values - {'?'}) == 1:
            column_types.append('invariant')
        else:
            column_types.append('variant')

    # Step 2: Identify contiguous invariant columns and their lengths
    current_invariant_block = None
    for col_idx in range(num_columns):
        if column_types[col_idx] == 'invariant':
            if current_invariant_block is None:
                current_invariant_block = {'start': col_idx, 'length': 1}
            else:
                current_invariant_block['length'] += 1
        else:
            if current_invariant_block:
                contiguous_invariant_blocks.append(current_invariant_block)
                current_invariant_block = None
    if current_invariant_block:
        contiguous_invariant_blocks.append(current_invariant_block)

    # Step 3: Log and optionally process blocks
    contiguous_invariant_blocks.sort(key=lambda x: x['length'], reverse=True)
    block_lengths = {}
    for block in contiguous_invariant_blocks:
        block_lengths.setdefault(block['length'], []).append(block)

    if partitioning_round == "max":
        add_breaks_terminal(alignment)  # Call the other function instead of inserting hashtags
    else:
        # Step 4: Track the positions for inserting hashtags
        hashtag_positions = []
        block_lengths_sorted = sorted(block_lengths.keys(), reverse=True)
        for block_length in block_lengths_sorted[:partitioning_round]:
            blocks = block_lengths[block_length]
            for block in blocks:
                start_idx = block['start']
                end_idx = start_idx + block['length'] - 1
                middle_idx = (start_idx + end_idx) // 2
                hashtag_positions.append(middle_idx)

        # Step 5: Insert hashtags
        for seq_id, seq in alignment.items():
            sorted_positions = sorted(hashtag_positions)
            shift = 0
            for middle_idx in sorted_positions:
                adjusted_idx = middle_idx + shift
                seq = seq[:adjusted_idx + 1] + '#' + seq[adjusted_idx + 1:]
                shift += 1
            alignment[seq_id] = seq

    # Step 6: Optionally log to CSV
    if log_csv_output:
        file_dir = os.path.dirname(csv_file_path)
        if not os.path.exists(file_dir) and file_dir:
            os.makedirs(file_dir)

        with open(csv_file_path, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Start', 'End', 'Length'])
            for block in contiguous_invariant_blocks:
                start_idx = block['start']
                end_idx = start_idx + block['length'] - 1
                writer.writerow([start_idx, end_idx, block['length']])
        print(f"Log of contiguous invariant blocks written to {csv_file_path}")

    return alignment, contiguous_invariant_blocks

def refinement_question2hyphen(alignment):
    """
    Replaces contiguous '?' characters flanked by '#' with '-' in the sequence alignment.
    This includes blocks surrounded by '#' as well as the first and last blocks that are flanked only on one side.

    Parameters:
        alignment (dict): Dictionary with sequence names as keys and sequences as values.

    Returns:
        dict: Modified alignment with '?' replaced by '-' where flanked by '#'.
    """
    # Iterate through each sequence in the alignment
    for name, sequence in alignment.items():
        # Replace '?' surrounded by '#' on both sides (internal blocks)
        modified_sequence = re.sub(r'(?<=#)\?+(?=#)', lambda m: '-' * len(m.group(0)), sequence)
        
        # Replace '?' in the first block (flanked by # on the right)
        modified_sequence = re.sub(r'^(\?+)(?=#)', lambda m: '-' * len(m.group(0)), modified_sequence)
        
        # Replace '?' in the last block (flanked by # on the left)
        modified_sequence = re.sub(r'(?<=#)(\?+)$', lambda m: '-' * len(m.group(0)), modified_sequence)

        # Update the alignment dictionary with the modified sequence
        alignment[name] = modified_sequence

    return alignment

def remove_columns_with_W(alignment: dict) -> dict:
    """
    Remove problematic columns in a DNA alignment:
    
    1. For each sequence, identify 15 contiguous 'W' or 'w' characters. 
       Mark those columns for deletion.
    2. Remove columns where all values are:
       - only 'W'
       - only 'W' and '?'
       - only 'W', '?', and '-'
    
    Args:
        alignment (dict): Dictionary of {sequence_id: sequence_string}

    Returns:
        dict: Cleaned alignment with columns removed
    """
    import numpy as np

    # Convert alignment to matrix
    seq_ids = list(alignment.keys())
    seqs = list(alignment.values())
    alignment_array = np.array([list(seq) for seq in seqs])
    n_rows, n_cols = alignment_array.shape

    columns_to_remove = set()

    # Step 1: Find 15 contiguous W/w in any sequence
    for row in alignment_array:
        upper_row = [char.upper() for char in row]
        for i in range(n_cols - 14):
            if all(base == 'W' for base in upper_row[i:i+15]):
                columns_to_remove.update(range(i, i + 15))

    # Step 2: Remove columns with only W / W+? / W+?+-
    for col_idx in range(n_cols):
        col_bases = set(alignment_array[:, col_idx].astype(str).flatten().tolist())
        col_bases_upper = {b.upper() for b in col_bases}
        if col_bases_upper.issubset({'W'}) or \
           col_bases_upper.issubset({'W', '?'}) or \
           col_bases_upper.issubset({'W', '?', '-'}):
            columns_to_remove.add(col_idx)

    # Remove columns
    columns_to_keep = sorted(set(range(n_cols)) - columns_to_remove)
    cleaned_array = alignment_array[:, columns_to_keep]

    # Convert back to dictionary
    cleaned_alignment = {
        seq_id: ''.join(cleaned_array[i]) for i, seq_id in enumerate(seq_ids)
    }

    return cleaned_alignment

def n2question_func(alignment: dict, leaves='all', log=False):
    """
    Replace all ambiguous nucleotides 'N' or 'n' with '?' in selected sequences.

    Parameters:
    - alignment (dict): Dictionary where keys are sequence names and values are DNA sequences (str).
    - leaves (str or list): Sequence name(s) to apply the replacement. Use 'all' to apply to all sequences.
    - log (bool): If True, also return a log of replaced blocks with their positions.

    Returns:
    - dict: Modified alignment with 'N'/'n' replaced with '?'.
    - list (optional): List of tuples (seq_name, start, end) for each replaced block.
    """
    if leaves == 'all':
        leaves_to_process = alignment.keys()
    elif isinstance(leaves, str):
        leaves_to_process = [leaves]
    else:
        leaves_to_process = leaves

    updated_alignment = {}
    replacement_log = []

    for name, seq in alignment.items():
        if name in leaves_to_process:
            new_seq = []
            i = 0
            while i < len(seq):
                if seq[i] in ('N', 'n'):
                    start = i
                    while i < len(seq) and seq[i] in ('N', 'n'):
                        i += 1
                    end = i - 1
                    replacement_log.append((name, start, end))
                    new_seq.extend(['?'] * (end - start + 1))
                else:
                    new_seq.append(seq[i])
                    i += 1
            updated_alignment[name] = ''.join(new_seq)
        else:
            updated_alignment[name] = seq

    if log:
        return updated_alignment, replacement_log
    return updated_alignment

##############################
# AUXILIARY FUNCTIONS TO LOG #
##############################

def compute_summary_after(alignment):
    num_seqs = len(alignment)

    # Transpose to columns
    columns = list(zip(*alignment.values()))
    
    # Count columns that contain pound signs
    total_pound = sum('#' in col for col in columns)
    
    # Alignment length excluding columns of only pound signs
    aln_length = sum(1 for col in columns if set(col) != {'#'})

    # Count nucleotide and gap characters
    total_nt = sum(c in "ACGTacgt" for seq in alignment.values() for c in seq)
    total_gaps = sum(seq.count("-") for seq in alignment.values())
    total_ns = sum(c in "Nn" for seq in alignment.values() for c in seq)
    total_qm = sum(seq.count("?") for seq in alignment.values())  # <- includes everything now

    # Count additional missing data as gap-only blocks between #
    missing_by_partition = 0
    for seq in alignment.values():
        parts = seq.split("#")
        for part in parts:
            if all(c == '-' for c in part):
                missing_by_partition += len(part)

    total_missing = total_qm + missing_by_partition

    return {
        "num_seqs": num_seqs,
        "aln_length": aln_length,
        "total_nt": total_nt,
        "total_gaps": total_gaps,
        "total_ns": total_ns,
        "total_qm": total_qm,
        "total_pound": total_pound,
        "missing_by_partition": missing_by_partition,
        "total_missing": total_missing
    }

def detect_fully_missing_partitions(alignment):
    """
    Logs all `?` and `-` from fully missing partitions.
    A fully missing partition is a region (between #) in which a sequence has only dashes.
    """
    log_entries = []
    total_question_marks = 0
    total_dash_from_missing_partitions = 0

    for seq_id, seq in alignment.items():
        total_question_marks += seq.count("?")

        parts = seq.split("#")
        col_index = 0  # absolute position tracker
        for i, part in enumerate(parts):
            if all(c == '-' for c in part):
                dash_count = len(part)
                total_dash_from_missing_partitions += dash_count
                start = col_index
                end = col_index + dash_count - 1
                log_entries.append(f"{seq_id}: partition {i} ({start}–{end}, length {dash_count}) fully missing (all '-')")
            col_index += len(part) + 1  # +1 for the '#' removed in split

    summary = (
        f"Total '?' characters: {total_question_marks}\n"
        f"Total '-' characters in fully missing partitions: {total_dash_from_missing_partitions}\n"
        f"Combined total: {total_question_marks + total_dash_from_missing_partitions}\n"
    )

    return summary + "\n" + "\n".join(log_entries)

#######################
# INTEGRATED FUNCTIONS #
#######################

def GB2MSA(input_file, 
           output_prefix, 
           delimiter=',', 
           write_names=True, 
           log=False, 
           orphan_threshold=6):
    """
    Complete GenBank-to-MSA pipeline:
    1. Downloads sequences from GenBank and aligns them by gene using MAFFT.
    2. Cleans the alignments by replacing internal missing data and removing empty columns.
    3. Applies GB2MSA_3 to replace selected gap and orphan nucleotide blocks with '?'.
    4. Applies GB2MSA_4 to replace blocks of 15 or more w/W (with or without interspersed gaps) with '?'.
    5. Replaces terminal '?' in sequences with '-'.
    6. Deletes intermediate files ending with '_aligned.fasta'.
    7. Optionally logs wall clock and CPU time to a log file named '<output_prefix>_log.txt'.
    """
    start_wall = time.time()
    start_cpu = time.process_time()

    # Step 1: Generate aligned FASTA files
    aligned_files = GB2MSA_1(input_file, output_prefix, delimiter=delimiter, write_names=write_names)

    # Step 2: Clean each aligned FASTA file
    cleaned_files = []
    for aligned_file in aligned_files:
        cleaned_file = GB2MSA_2(aligned_file)
        cleaned_files.append(cleaned_file)

    # Step 3: Apply GB2MSA_3 to clean orphan gap/nucleotide blocks
    for cleaned_file in cleaned_files:
        records = list(SeqIO.parse(cleaned_file, "fasta"))
        alignment_dict = {record.id: str(record.seq) for record in records}
        updated_dict = GB2MSA_3(alignment_dict, orphan_threshold=orphan_threshold, log=log)
        
        # Step 4: Apply GB2MSA_4 to handle w/W blocks
        updated_dict = GB2MSA_4(updated_dict)

        updated_records = []
        for record in records:
            record.seq = Seq(updated_dict[record.id])
            updated_records.append(record)
        with open(cleaned_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "fasta")

    # Step 5: Replace terminal '?' with '-' in each sequence
    for cleaned_file in cleaned_files:
        records = list(SeqIO.parse(cleaned_file, "fasta"))
        updated_records = []
        for record in records:
            seq = str(record.seq)
            left = len(seq) - len(seq.lstrip('?'))
            right = len(seq) - len(seq.rstrip('?'))
            new_seq = '-' * left + seq[left:len(seq)-right] + '-' * right if right > 0 else '-' * left + seq[left:]
            record.seq = Seq(new_seq)
            updated_records.append(record)
        with open(cleaned_file, "w") as out_handle:
            SeqIO.write(updated_records, out_handle, "fasta")

    # Step 6: Delete intermediate *_aligned.fasta files
    for aligned_file in aligned_files:
        if aligned_file.endswith("_aligned.fasta") and os.path.exists(aligned_file):
            os.remove(aligned_file)

    # Step 7: Log timing if requested
    if log:
        end_wall = time.time()
        end_cpu = time.process_time()
        wall_time = end_wall - start_wall
        cpu_time = end_cpu - start_cpu

        log_file = f"{output_prefix}_log.txt"
        with open(log_file, "a") as lf:
            lf.write(f"--- GB2MSA run for '{output_prefix}' ---\n")
            lf.write(f"Wall clock time: {wall_time:.2f} seconds\n")
            lf.write(f"CPU time: {cpu_time:.2f} seconds\n\n")

    return cleaned_files

def addSeq(
    alignment,
    new_seqs,
    output,
    write_names=True,
    orphan_threshold=0,
    log=False,
    n2question=None,
    gaps2question=None
):
    """
    Add new sequences to an existing alignment using MAFFT, clean and standardize the result.

    Steps performed:
    1. Remove '#' columns from original alignment.
    2. Align new sequences with MAFFT using --add.
    3. Trim orphan nucleotide blocks from new sequences.
    4. Remove short DNA blocks near # in first and last partitions.
    5. Replace terminal '-' with '?'.
    6. Optionally replace all N/n with '?' in selected sequences.
    7. Optionally replace long gap blocks with '?'.
    8. Reinsert '#' columns.
    9. Replace all-'?' blocks between '#' with '-'.
    10. Write the result and an optional log file.

    Parameters:
        alignment (str or dict): Existing alignment file path (FASTA) or dictionary {id: sequence}.
        new_seqs (str or dict): New sequences to add (FASTA path or dict {id: sequence}).
        output (str): Path to write the updated alignment in FASTA format.
        write_names (bool): Whether to write a _terminal_names.txt file listing sequence IDs.
        orphan_threshold (int): Threshold to detect and remove orphan DNA blocks.
        log (bool): If True, writes a log file with trimming and runtime information.
        n2question (str, list or None): Replace 'N/n' with '?' in specific sequences:
            - 'all': apply to all sequences
            - str: apply to a single sequence ID
            - list: apply to listed sequence IDs
        gaps2question (int or None): Replace contiguous gap blocks larger than this threshold with '?'. Only applied to added sequences.

    Returns:
        None

    Example usage:
        addSeq("alignment.fasta", "new.fasta", "updated.fasta", n2question="seq123", log=True)
        addSeq(alignment_dict, new_dict, "out.fas", n2question='all')
    """

    # Start timing the execution
    start_time = time.time()
    temp_files_to_remove = []  # Temporary files to be removed after execution
    log_lines = [] 

    # Log the function call and parameters for reproducibility
    if log:
        cmd_used = f"addSeq(alignment=..., new_seqs=..., output='{output}', write_names={write_names}, orphan_threshold={orphan_threshold}, log={log}, n2question={n2question}, gaps2question={gaps2question})"
        log_lines.append(f"Command used: {cmd_used}")
        log_lines.append("")

    # === Step 1: Load and clean the alignment ===
    def write_dict_to_temp_fasta(seq_dict):
        records = [SeqRecord(Seq(seq), id=str(seq_id), description="") for seq_id, seq in seq_dict.items()]
        tmp = tempfile.NamedTemporaryFile("w+", delete=False)
        SeqIO.write(records, tmp, "fasta")
        tmp.close()
        return tmp.name

    if isinstance(alignment, dict):
        alignment_path = write_dict_to_temp_fasta(alignment)
        temp_files_to_remove.append(alignment_path)
    elif isinstance(alignment, str):
        alignment_path = alignment
    else:
        raise ValueError("alignment must be a FASTA file path or a dictionary")

    records = list(SeqIO.parse(alignment_path, "fasta"))
    if not records:
        raise ValueError("Input alignment is empty or not found")

    aln_len = len(records[0].seq)
    # Identify columns that are '#' characters to temporarily remove them for alignment
    pound_cols = [i for i in range(aln_len) if any(rec.seq[i] == '#' for rec in records)]

    # Log input alignment info
    if log:
        log_lines.append(f"Input alignment: {len(records)} sequences")
        log_lines.append(f"Input alignment: {len(pound_cols)} # columns")

    def remove_cols(seq, cols):
        return ''.join(seq[i] for i in range(len(seq)) if i not in cols)

    # Remove '#' columns
    aln_no_pound = [
        SeqRecord(Seq(remove_cols(str(rec.seq), pound_cols)), id=rec.id, description="")
        for rec in records
    ]
    # Write cleaned alignment to a temporary file
    with tempfile.NamedTemporaryFile("w+", delete=False) as aln_tmp:
        SeqIO.write(aln_no_pound, aln_tmp, "fasta")
        aln_path = aln_tmp.name
        temp_files_to_remove.append(aln_path)

    # === Step 2: Load new sequences ===
    if isinstance(new_seqs, dict):
        new_seqs_path = write_dict_to_temp_fasta(new_seqs)
        new_seq_count = len(new_seqs)
        temp_files_to_remove.append(new_seqs_path)
    elif isinstance(new_seqs, str):
        new_seq_count = sum(1 for _ in SeqIO.parse(new_seqs, "fasta"))
        new_seqs_path = new_seqs
    else:
        raise ValueError("new_seqs must be a FASTA file path or a dictionary")
    # Log new sequence info
    if log:
        log_lines.append(f"Input new sequences: {new_seq_count} sequences")

    # === Step 3: Align new sequences with MAFFT ===
    with tempfile.NamedTemporaryFile("w+", delete=False) as out_tmp:
        out_path = out_tmp.name
        temp_files_to_remove.append(out_path)

    try:
        subprocess.run(
            ['mafft', '--add', new_seqs_path, '--keeplength', '--preservecase', aln_path],
            check=True,
            stdout=open(out_path, 'w'),
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        for file in temp_files_to_remove:
            try:
                os.remove(file)
            except Exception:
                pass
        raise RuntimeError(f"MAFFT failed!\nCommand: {e.cmd}\nExit status: {e.returncode}\nMAFFT error output:\n{e.stderr}")

    mafft_aligned_records = list(SeqIO.parse(out_path, "fasta"))
    original_ids = {rec.id for rec in records}
    new_records = [rec for rec in mafft_aligned_records if rec.id not in original_ids]


    def replace_gap_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        replaced_log = []
        i = 0
        while i < len(seq_list):
            if seq_list[i] == '-':
                start = i
                while i < len(seq_list) and seq_list[i] == '-':
                    i += 1
                if (i - start) > threshold:
                    for j in range(start, i):
                        seq_list[j] = '?'
                    if seq_id:
                        replaced_log.append(f"{seq_id}: {i - start} contiguous '-' replaced with '?' at {start}–{i}")
            else:
                i += 1
        return ''.join(seq_list), replaced_log

    def find_dna_blocks(seq_list, start, end):
        blocks = []
        i = start
        while i < end:
            if seq_list[i] not in "-?#":
                s = i
                while i < end and seq_list[i] not in "-?#":
                    i += 1
                e = i
                blocks.append((s, e))
            else:
                i += 1
        return blocks

    def trim_orphan_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        trimmed_log = []

        def find_blocks(seq_list):
            blocks = []
            i = 0
            while i < len(seq_list):
                if seq_list[i] not in "-?#":
                    start = i
                    while i < len(seq_list) and seq_list[i] not in "-?#":
                        i += 1
                    end = i
                    blocks.append((start, end))
                else:
                    i += 1
            return blocks

        changed = True
        while changed:
            changed = False

            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                first_start, first_end = blocks[0]
                next_start = blocks[1][0]
                gap_count = seq_list[first_end:next_start].count('-') + seq_list[first_end:next_start].count('?')
                size = first_end - first_start
                if size < threshold and gap_count > threshold:
                    deleted = ''.join(seq_list[first_start:first_end])
                    seq_list[first_start:first_end] = ['-'] * size
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Left {first_start}–{first_end} (size={size}, '{deleted}')")
                    changed = True
                    continue
                break

            while True:
                blocks = find_blocks(seq_list)
                if len(blocks) < 2:
                    break
                last_start, last_end = blocks[-1]
                prev_end = blocks[-2][1]
                gap_count = seq_list[prev_end:last_start].count('-') + seq_list[prev_end:last_start].count('?')
                size = last_end - last_start
                if size < threshold and gap_count > threshold:
                    deleted = ''.join(seq_list[last_start:last_end])
                    seq_list[last_start:last_end] = ['-'] * size
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Right {last_start}–{last_end} (size={size}, '{deleted}')")
                    changed = True
                    continue
                break

        if pound_cols:
            first_hash = min(pound_cols)
            last_hash = max(pound_cols)
            blocks = find_dna_blocks(seq_list, 0, first_hash)
            if len(blocks) == 1:
                s, e = blocks[0]
                if e == first_hash and (e - s) < threshold:
                    deleted = ''.join(seq_list[s:e])
                    seq_list[s:e] = ['-'] * (e - s)
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Left {s}–{e} (size={e - s}, '{deleted}')")
            blocks = find_dna_blocks(seq_list, last_hash + 1, len(seq_list))
            if len(blocks) == 1:
                s, e = blocks[0]
                if s == last_hash + 1 and (e - s) < threshold:
                    deleted = ''.join(seq_list[s:e])
                    seq_list[s:e] = ['-'] * (e - s)
                    if seq_id:
                        trimmed_log.append(f"{seq_id}: Right {s}–{e} (size={e - s}, '{deleted}')")

        return ''.join(seq_list), trimmed_log

    trimmed_new_records = []
    all_trim_logs = []
    gaps2q_log = []  # This will no longer be used for logging replaced gaps

    for rec in new_records:
        trimmed_seq, seq_log = trim_orphan_blocks(str(rec.seq), orphan_threshold, seq_id=rec.id)
        # Remove gaps2question here to not log replaced gaps multiple times
        trimmed_new_records.append(SeqRecord(Seq(trimmed_seq), id=rec.id, description=""))
        all_trim_logs.extend(seq_log)

    processed_records = [rec for rec in mafft_aligned_records if rec.id in original_ids] + trimmed_new_records

    updated_records = []
    for rec in processed_records:
        seq_chars = list(str(rec.seq))
        for i in range(len(seq_chars)):
            if seq_chars[i] == '-':
                seq_chars[i] = '?'
            else:
                break
        for i in range(len(seq_chars) - 1, -1, -1):
            if seq_chars[i] == '-':
                seq_chars[i] = '?'
            else:
                break
        updated_records.append(SeqRecord(Seq(''.join(seq_chars)), id=rec.id, description=""))

    n2q_log = []
    if n2question:
        if isinstance(n2question, str) and n2question != 'all':
            target_ids = {n2question}
        elif isinstance(n2question, list):
            target_ids = set(n2question)
        elif n2question == 'all':
            target_ids = {rec.id for rec in updated_records}
        else:
            target_ids = set()

        for rec in updated_records:
            if rec.id in target_ids:
                seq_str = str(rec.seq)
                count_n = seq_str.count('N') + seq_str.count('n')
                if count_n > 0:
                    rec.seq = Seq(seq_str.replace('N', '?').replace('n', '?'))
                    n2q_log.append(f"{rec.id}: {count_n} N/n replaced with ?")

    # === Now apply gaps2question as the very last step on added sequences ===
    def replace_gap_blocks(seq, threshold, seq_id=None):
        seq_list = list(seq)
        replaced_log = []
        i = 0
        while i < len(seq_list):
            if seq_list[i] == '-':
                start = i
                while i < len(seq_list) and seq_list[i] == '-':
                    i += 1
                if (i - start) > threshold:
                    for j in range(start, i):
                        seq_list[j] = '?'
                    if seq_id:
                        replaced_log.append(f"{seq_id}: {i - start} contiguous '-' replaced with '?' at {start}–{i}")
            else:
                i += 1
        return ''.join(seq_list), replaced_log

    # Apply gaps2question only if specified
    if gaps2question:
        updated_records_dict = {rec.id: rec for rec in updated_records}
        for rec in trimmed_new_records:
            seq = str(updated_records_dict[rec.id].seq)
            new_seq, _ = replace_gap_blocks(seq, gaps2question, seq_id=rec.id)
            updated_records_dict[rec.id].seq = Seq(new_seq)
        updated_records = list(updated_records_dict.values())

    final_records = []
    for rec in updated_records:
        seq_list = list(str(rec.seq))
        for col in sorted(pound_cols):
            seq_list.insert(col, '#')
        final_records.append(SeqRecord(Seq(''.join(seq_list)), id=rec.id, description=""))

    def process_blocks(seq):
        seq_chars = list(seq)
        blocks = []
        start = 0
        for i, c in enumerate(seq_chars):
            if c == '#':
                blocks.append((start, i))
                start = i + 1
        blocks.append((start, len(seq_chars)))
        for (start, end) in blocks:
            if all(seq_chars[i] == '?' for i in range(start, end)):
                for i in range(start, end):
                    seq_chars[i] = '-'
        return ''.join(seq_chars)

    final_output = [
        SeqRecord(Seq(process_blocks(str(rec.seq))), id=rec.id, description="")
        for rec in final_records
    ]

    # --- New function to find contiguous '?' blocks ---
    def find_question_blocks(seq):
        blocks = []
        seq_len = len(seq)
        i = 0
        while i < seq_len:
            if seq[i] == '?':
                start = i
                while i < seq_len and seq[i] == '?':
                    i += 1
                end = i
                blocks.append((start, end))
            else:
                i += 1
        return blocks

    # Identify gap ('?') blocks in added sequences for final logging
    gap_question_blocks_log = []
    added_ids = {rec.id for rec in trimmed_new_records}  # IDs of added sequences

    for rec in final_output:
        if rec.id in added_ids:
            q_blocks = find_question_blocks(str(rec.seq))
            for (start, end) in q_blocks:
                length = end - start
                gap_question_blocks_log.append(f"{rec.id}: ? block at positions {start}-{end} (length={length})")

    SeqIO.write(final_output, output, "fasta")
    if write_names:
        with open(output + "_terminal_names.txt", "w") as f:
            for rec in final_output:
                f.write(rec.id + "\n")

    if log:
        elapsed = time.time() - start_time
        log_lines.append(f"Final alignment: {len(final_output)} sequences")
        log_lines.append(f"Final alignment: {len(final_output[0].seq)} columns")
        log_lines.append(f"Final alignment: {sum(1 for i in range(len(final_output[0].seq)) if any(rec.seq[i] == '#' for rec in final_output))} # columns")
        log_lines.append("")
        log_lines.append("Trimmed orphan blocks from new sequences:")
        log_lines.extend(all_trim_logs or ["None"])
        if gap_question_blocks_log:
            log_lines.append("")
            log_lines.append("Gap block replacements:")
            log_lines.extend(gap_question_blocks_log)
        if n2q_log:
            log_lines.append("")
            log_lines.append("N/n to ? replacements:")
            log_lines.extend(n2q_log)
        log_lines.append("")
        log_lines.append(f"Runtime: {elapsed:.2f} seconds")
        with open(output + ".log", "w") as log_file:
            log_file.write("\n".join(log_lines))

    for file in temp_files_to_remove:
        try:
            os.remove(file)
        except Exception:
            pass

def prepDyn(input_file=None, 
            GB_input=None,
            input_format="fasta",
            MSA=False,
            output_file=None,
            output_format="fasta",
            log=False,
            # Trimming parameters
            orphan_method=None,
            orphan_threshold=10,
            percentile=25,
            del_inv=True,
            # Missing data parameters
            internal_method=None,
            internal_column_ranges=None,
            internal_leaves="all",
            internal_threshold=None,
            n2question=None,
            # Partitioning parameters
            partitioning_round=0
            ):
    """
    Preprocess missing data for dynamic homology in PhyG. First, columns containing 
    only gaps, orphan nucleotides, and invariant columns can be trimmed. Second, 
    missing data is coded with question marks. Third, partitions are delimited in 
    highly conserved regions.

    Args:
        input_file (str): Path to the input alignment file or directory. Ignored if GB_input is provided.
        GB_input (str): Path to a CSV/TSV file containing GenBank accession numbers. If provided,
                        sequences will be downloaded from GenBank and aligned before preprocessing.
        input_format (str): Format of the input alignment. Options: 'fasta' (default), 
                            'clustal', 'phylip', or any format accepted by Biopython. 
        MSA (bool): Whether to perform MSA if input sequences specified in input_file are unaligned
        orphan_method (str): The trimming method. By default, trimming orphan nucleotides
                             is not performed. Options:
                            - 'auto': trim using the 25th percentile;
                            - 'semi': trim with a manual threshold.
        orphan_threshold (int): Threshold used to trim orphan nucleotides if orphan_method = 'semi'.
        percentile (float): Used with orphan_method = 'auto' to define trimming threshold.
        del_inv (bool): Whether to trim invariant terminal columns. Default is True.
        internal_method (str): Defines how to identify internal missing data. Automatic identificaton
                               of missing data is made if GB_input is provided. Otherwise, naive 
                               options to identify internal missing data are:
                               - "manual": Use column ranges;
                               - "semi": Use a threshold for gaps.
        internal_column_ranges (list): Column ranges (inclusive) if internal_method = 'manual'.
        internal_leaves (str or list): Sequences to apply internal missing data replacement 
                                       if internal_method is not "None".
        internal_threshold (int): Used with internal_method = 'semi' to define gap threshold.
                                  Contiguous '-' larger than the threshold are replaced with '?'.
        partitioning_round (int): Number of partitioning round. Invariant regions are sorted by length
                                  in descendant order and the n-largest block(s) partitioned using '#'.
                                  If "max" is specified, pound signs are inserted arund all blocks of 
                                  missing data. 
        output_file (str): Custom prefix for output files. If None, base_name from input_file is used.
        output_format (str): Output format. Default is 'fasta'.
        log (bool): Whether to write a log with wall-clock time. Default is False.
        n2question (str or list): If specified, replaces ambiguous nucleotide 'N' or 'n' with '?'. If None (default), n2question is not performed. If 'all', apply to all sequences. If you want to apply to only one sequence, write the name of this sequence. If you want to apply to multiple sequences but no all, wrie the list of sequences.
                                  
    Returns:
        dict: The preprocessed unaligned sequences.
    """

    # Start timers if logging is enabled
    if log:
        start_wall_time = time.time()
        start_cpu_time = time.process_time()

    # Step 1: Run GB2MSA if GenBank input is provided
    if GB_input is not None:
        print("Running GB2MSA on GenBank input...")
        gb_output_prefix = output_file if output_file else "output"
        cleaned_files = GB2MSA(GB_input, output_prefix=gb_output_prefix, log=False)

        for file in cleaned_files:
            gene_name = os.path.splitext(os.path.basename(file))[0].replace("output_", "")
            if output_file:
                specific_output_prefix = f"{output_file}_{gene_name}"
            else:
                specific_output_prefix = f"{gene_name}"

            alignment = AlignIO.read(file, "fasta")
            alignment_dict = {record.id: str(record.seq) for record in alignment}
            
            prepDyn(input_file=alignment_dict,
                    input_format="dict",
                    orphan_method=orphan_method,
                    orphan_threshold=orphan_threshold,
                    percentile=percentile,
                    del_inv=del_inv,
                    internal_method=internal_method,
                    internal_column_ranges=internal_column_ranges,
                    internal_leaves=internal_leaves,
                    internal_threshold=internal_threshold,
                    n2question=n2question,
                    partitioning_round=partitioning_round,
                    output_format=output_format,
                    log=log,
                    output_file=specific_output_prefix)

        return

    # Step 2: If a folder is provided, process each alignment inside
    if isinstance(input_file, str) and os.path.isdir(input_file):
        for file_name in os.listdir(input_file):
            if file_name.endswith(f".{input_format}"):
                file_path = os.path.join(input_file, file_name)
                prepDyn(file_path,
                      input_format=input_format,
                      MSA=MSA,
                      orphan_method=orphan_method,
                      orphan_threshold=orphan_threshold,
                      percentile=percentile,
                      del_inv=del_inv,
                      internal_method=internal_method,
                      internal_column_ranges=internal_column_ranges,
                      internal_leaves=internal_leaves,
                      internal_threshold=internal_threshold,
                      n2question=n2question,
                      output_format=output_format)
        return

    # Step 3: Read and process alignment
    if isinstance(input_file, dict):
        alignment = input_file
    else:
        alignment = AlignIO.read(input_file, input_format)
        alignment = {record.id: str(record.seq) for record in alignment}

    # Optional MAFFT alignment
    if MSA:
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_in:
            for k, v in alignment.items():
                tmp_in.write(f">{k}\n{v}\n")
            tmp_in_path = tmp_in.name

        tmp_out_path = tmp_in_path + "_aligned.fasta"
        try:
            subprocess.run(["mafft", "--auto", tmp_in_path], stdout=open(tmp_out_path, "w"), stderr=subprocess.DEVNULL, check=True)
        except subprocess.CalledProcessError:
            raise RuntimeError("MAFFT alignment failed.")

        alignment = AlignIO.read(tmp_out_path, "fasta")
        alignment = {record.id: str(record.seq) for record in alignment}

        os.remove(tmp_in_path)
        os.remove(tmp_out_path)

    # Summary of characteristics before preprocessing
    if log:
        num_seqs = len(alignment)
        aln_length = len(next(iter(alignment.values())))
        total_nt = sum(c.upper() in "ACGT" for seq in alignment.values() for c in seq)
        total_gaps = sum(seq.count("-") for seq in alignment.values())
        total_ns = sum(c in "Nn" for seq in alignment.values() for c in seq)

    # 3.1 Remove columns with gaps in all leaves
    remove_all_gap_columns(alignment)
    
    # 3.2 Trim orphan nucleotides
    orphan_log = None
    if orphan_method == "percentile":
        orphan_threshold = calculate_orphan_threshold_from_percentile(alignment, percentile, terminal_only=True)
        if log:
            alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=True)
        else:
            alignment = delete_orphan_nucleotides2(alignment, orphan_threshold)
    elif orphan_method == "semi":
        if log:
            alignment, orphan_log = delete_orphan_nucleotides2(alignment, orphan_threshold, log_changes=True)
        else:
            alignment = delete_orphan_nucleotides2(alignment, orphan_threshold)


    # 3.3 Replace terminal gaps with ?
    alignment = replace_terminal_gaps_dict(alignment)

    # 3.4 Trim invariant columns
    removed_cols = []
    if del_inv:
        alignment = remove_non_informative_positions(alignment, removed_indices=removed_cols)

    alignment = replace_terminal_gaps_dict(alignment)

    # 3.5 Replace internal gaps with ?
    if internal_method == "manual":
        alignment = replace_dashes_with_question_marks(alignment=alignment, 
                                                       internal_column_ranges=internal_column_ranges, 
                                                       internal_leaves=internal_leaves, 
                                                       internal_method="manual")
    elif internal_method == "semi":
        alignment = replace_dashes_with_question_marks(alignment=alignment, 
                                                       internal_leaves=internal_leaves, 
                                                       internal_method="semi", 
                                                       internal_threshold=internal_threshold)

    # 3.6 Replace ambiguous nucleotides N/n with ?
    n_blocks = []
    if n2question is not None:
        alignment, n_blocks = n2question_func(alignment, leaves=n2question, log=True)


    # 3.7 Successive partition
    classify_and_insert_hashtags(alignment, 
                                 partitioning_round=partitioning_round)
    refinement_question2hyphen(alignment)
    alignment = remove_columns_with_W(alignment)

    # Step 4: Write output file
    records = [SeqRecord(Seq(seq), id=key, description="") for key, seq in alignment.items()]
    base_name = os.path.splitext(os.path.basename(input_file))[0] if isinstance(input_file, str) else "alignment"

    if output_file:
        output_prefix = f"{output_file}"
    else:
        output_prefix = base_name

    output_path = f"{output_prefix}_preprocessed.{output_format}"
    with open(output_path, "w") as output_handle:
        SeqIO.write(records, output_handle, output_format)

    # Step 5: Write log
    if log:
        end_wall_time = time.time()
        end_cpu_time = time.process_time()
        wall_time = end_wall_time - start_wall_time
        cpu_time = end_cpu_time - start_cpu_time
        
        # Build the command string with all parameters and their current values
        cmd_parts = ["prepDyn("]
        params = {
            "input_file": input_file,
            "GB_input": GB_input,
            "input_format": input_format,
            "MSA": MSA,
            "output_file": output_file,
            "output_format": output_format,
            "log": log,
            "orphan_method": orphan_method,
            "orphan_threshold": orphan_threshold,
            "percentile": percentile,
            "del_inv": del_inv,
            "internal_method": internal_method,
            "internal_column_ranges": internal_column_ranges,
            "internal_leaves": internal_leaves,
            "internal_threshold": internal_threshold,
            "n2question": n2question,
            "partitioning_round": partitioning_round
        }
        param_strs = []
        for k, v in params.items():
            if v is None:
                param_strs.append(f"{k}=None")
            elif isinstance(v, str):
                param_strs.append(f"{k}='{v}'")
            else:
                param_strs.append(f"{k}={repr(v)}")
        cmd_parts.append(", ".join(param_strs))
        cmd_parts.append(")")
        cmd_line = "".join(cmd_parts)
        
        log_path = f"{output_prefix}_log.txt"
        with open(log_path, "w") as log_file:
            log_file.write("--- Command used ---\n")
            log_file.write(f"{cmd_line}\n\n")

            log_file.write("--- Step 1: Summary before preprocessing ---\n")
            log_file.write(f"No. sequences: {num_seqs}\n")
            log_file.write(f"No. columns: {aln_length}\n")
            log_file.write(f"Total no. nucleotides (A/C/G/T only): {total_nt} bp\n")
            log_file.write(f"Total no. gaps (-): {total_gaps}\n")
            log_file.write(f"Total no. IUPAC N: {total_ns}\n\n")

            # log trimming: invariants
            if del_inv:
                log_file.write("--- Step 2: Trimming (invariant columns) ---\n")
                log_file.write(f"{removed_cols}\n\n")
            # log trimming: orphans
            if orphan_log:
                log_file.write("--- Step 2: Trimming (orphan nucleotides) ---\n")
                log_file.write(f"{orphan_log}\n\n")
            
            # log missing: n2question
            if n_blocks:
                log_file.write("--- Step 3: Missing data identification (Ns replaced with '?') ---\n")
                for seq_name, start, end in n_blocks:
                    log_file.write(f"{seq_name}: {start}–{end}\n")
                log_file.write("\n")
            
            # log terminal ? detection
            missing_partition_log = detect_fully_missing_partitions(alignment)
            if missing_partition_log:
                log_file.write("--- Step 3: Missing data identification ---\n")
                log_file.write(f"{missing_partition_log}\n\n")


            # log partitioning (#)
            if partitioning_round > 0:
                log_file.write("--- Step 4: Partitioning (columns with # inserted) ---\n")
                # Transpose alignment to columns
                columns = list(zip(*alignment.values()))
                pound_indices = [i for i, col in enumerate(columns) if '#' in col]
                log_file.write(f"{pound_indices}\n\n")
            
            # log preprocessed summary            
            summary_post = compute_summary_after(alignment)
            log_file.write("--- Summary after preprocessing ---\n")
            log_file.write(f"No. sequences: {summary_post['num_seqs']}\n")
            log_file.write(f"No. columns: {summary_post['aln_length']}\n")
            log_file.write(f"No. pound sign columns (#): {summary_post['total_pound']}\n")
            log_file.write(f"Total no. nucleotides (A/C/G/T): {summary_post['total_nt']} bp\n")
            log_file.write(f"Total no. gaps (-): {summary_post['total_gaps']}\n")
            log_file.write(f"Total no. IUPAC N: {summary_post['total_ns']}\n")
            log_file.write(f"Total no. missing values (?): {summary_post['total_missing']}\n\n")
            
            # log run time
            log_file.write("--- Run time ---\n")
            log_file.write(f"Wall-clock time: {wall_time:.8f} seconds\n")
            log_file.write(f"CPU time: {cpu_time:.8f} seconds\n")

    return alignment
        