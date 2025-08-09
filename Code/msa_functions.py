# ===============================================================
# msa_functions.py: Helper functions for Vibrio genome multiple sequence alignment (msa) analysis
# ===============================================================
# This file contains all of the functions for use in the main coding file,
# "vibrio-analysis.ipynb". These functions perform sequence parsing,
# alignment statistics, k-mer scanning, variable region extraction,
# GC content calculations, and metadata parsing for downstream analyses.
# ===============================================================

# Import required Python packages
# ----------------------------
import os
import csv
import numpy as np
from collections import Counter
from Bio import SeqIO          # Biopython's FASTA/sequence parser
from glob import glob
from itertools import combinations
import re
import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO        # For reading alignment files
from Bio.Align import MultipleSeqAlignment
import math                    # Required for shannon_entropy()

# Function 1: count_zeroes_in_consensus(block)
# -------------------------------------------------
# Builds a consensus matrix for a set of sequences (block) and counts the
# number of bases (A, T, C, G, or 'other') that have zero counts in each column.
# Used as a measure of variability or gaps across positions.
def count_zeroes_in_consensus(block):
    block_size = len(block[0])
    base_order = ['A', 'T', 'C', 'G']
    num_zeroes = 0

    for i in range(block_size):
        counts = Counter(seq[i].upper() for seq in block)
        # Count missing occurrences of standard bases
        for base in base_order:
            if counts.get(base, 0) == 0:
                num_zeroes += 1
        # Count "other" characters with zero count
        others = set(counts) - set(base_order)
        if sum(counts[b] for b in others) == 0:
            num_zeroes += 1

    return num_zeroes

# Function 2: is_perfect_match(block)
# -------------------------------------------------
# Returns True if all sequences in a block are identical.
def is_perfect_match(block):
    return all(seq == block[0] for seq in block)

# Function 3: calculate_pairwise_identity(sequences)
# -------------------------------------------------
# Computes min, max, and average pairwise identity (%) between sequences.
# Ignores positions containing gaps ("-") or ambiguous bases ("N").
def calculate_pairwise_identity(sequences):
    if len(sequences) < 2:
        return 100.0, 100.0, 100.0  # If only one sequence, identity is 100%

    # Ensure equal lengths for all sequences
    min_len = min(len(seq) for seq in sequences)
    sequences = [seq[:min_len] for seq in sequences]
    identities = []

    # Compare each pair of sequences
    for seq1, seq2 in combinations(sequences, 2):
        matches = 0
        valid_positions = 0
        for nt1, nt2 in zip(seq1, seq2):
            if nt1 in ['-', 'N'] or nt2 in ['-', 'N']:
                continue
            valid_positions += 1
            if nt1 == nt2:
                matches += 1
        if valid_positions > 0:
            identity = (matches / valid_positions) * 100
            identities.append(identity)

    # Calculate summary stats
    if identities:
        avg_identity = sum(identities) / len(identities)
        min_identity = min(identities)
        max_identity = max(identities)
    else:
        avg_identity = min_identity = max_identity = 0.0

    return min_identity, max_identity, avg_identity

# Function 4: process_alignment(aligned_path, block_size=20)
# -------------------------------------------------
# Processes all *_aligned.fasta files in a directory.
# Calculates various sequence conservation metrics and returns results per file.
def process_alignment(aligned_path, block_size=20):
    fasta_files = glob(os.path.join(aligned_path, "*_aligned.fasta"))
    all_results = []

    for file_path in fasta_files:
        sequences = [str(record.seq).upper() for record in SeqIO.parse(file_path, "fasta")]
        if not sequences:
            continue

        # Align to shortest sequence length
        min_len = min(len(seq) for seq in sequences)
        sequences = [seq[:min_len] for seq in sequences]

        zero_counts = []
        perfect_match_indices = []

        # Sliding window analysis
        for i in range(min_len - block_size + 1):
            block = [seq[i:i+block_size] for seq in sequences]
            zero_count = count_zeroes_in_consensus(block)
            zero_counts.append(zero_count)
            if is_perfect_match(block):
                perfect_match_indices.append(i)

        # Compute pairwise identity stats
        min_id, max_id, avg_id = calculate_pairwise_identity(sequences)

        if zero_counts:
            zero_counts = np.array(zero_counts)
            max_val = np.max(zero_counts)
            mean_val = np.mean(zero_counts)
            var_val = np.var(zero_counts)
            count_76 = np.sum(zero_counts == 76)
            num_perfect = len(perfect_match_indices)
            max_distance = np.max(np.diff(perfect_match_indices)) if num_perfect >= 2 else 0
            num_windows = min_len - block_size + 1
            norm_perfect = num_perfect / num_windows if num_windows > 0 else 0

            result = {
            # The name of the aligned FASTA file being analyzed, e.g. "GC_00000990_aligned.fasta".
            'filename': os.path.basename(file_path),
                
            #The minimum sequence length among the sequences in the alignment.
            'alignment_length': min_len,

            # The maximum number of zeroes in any consensus column in the alignment.
            'max_zeroes': int(max_val),

            # The average number of zeroes per consensus matrix column (across the alignment).
            'mean_zeroes': float(mean_val),

            # The variance of zero counts across the matrix.
            'variance_zeroes': float(var_val),

            # The number of contiguous blocks of length 76 base pairs with all-zero consensus values.
            'count_76_zero_blocks': int(count_76),

            # The total number of perfectly conserved 20-mers (regions with no sequence variation across all strains).
            'num_perfect_matches': int(num_perfect),

            # The longest distance (in bp) between any two conserved 20-mer regions.
            'max_dist_between_matches': int(max_distance),

            # This is a normalized count of perfect 20-mers, scaled by alignment length. 
            'normalized_perfect_match_rate': round(norm_perfect, 4),

            #  Lowest similarity between any two strains.
            'min_percent_identity': round(min_id, 2),

            # Highest similarity between any two strains.
            'max_percent_identity': round(max_id, 2),

            # Mean similarity across all pairwise comparisons.
            'avg_percent_identity': round(avg_id, 2)
            }

            all_results.append(result)

    return all_results

# Function 5: save_to_csv(results, output_file)
# -------------------------------------------------
# Saves a list of dictionaries (alignment results) to a CSV file.
def save_to_csv(results, output_file="alignment_stats.csv"):
    if not results:
        print("No results to save.")
        return
    fieldnames = list(results[0].keys())
    with open(output_file, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"Results saved to {output_file}")

# Function 6: get_perfect_match_indices(sequences, window_size)
# -------------------------------------------------
# Finds start positions of perfectly conserved k-mers (default = 20 bp).
def get_perfect_match_indices(sequences, window_size=20):
    indices = []
    for i in range(len(sequences[0]) - window_size + 1):
        block = [seq[i:i+window_size] for seq in sequences]
        if all(seq == block[0] for seq in block):
            indices.append(i)
    return indices

# Function 7: find_variable_regions(...)
# -------------------------------------------------
# Given perfect match positions, finds variable regions between them.
# Only keeps regions within the specified length bounds.
def find_variable_regions(perfect_match_indices, alignment_length, kmer_size=20,
                          min_length=200, max_length=400):
    variable_regions = []

    def keep(start, end):
        length = end - start + 1
        return min_length <= length <= max_length

    if not perfect_match_indices:
        if keep(0, alignment_length - 1):
            variable_regions.append((0, alignment_length - 1))
        return variable_regions

    # Check start region
    if perfect_match_indices[0] > 0:
        start, end = 0, perfect_match_indices[0] - 1
        if keep(start, end):
            variable_regions.append((start, end))

    # Check between perfect match regions
    for i in range(len(perfect_match_indices) - 1):
        start = perfect_match_indices[i] + kmer_size
        end = perfect_match_indices[i + 1] - 1
        if keep(start, end):
            variable_regions.append((start, end))

    # Check end region
    last_end = perfect_match_indices[-1] + kmer_size - 1
    if last_end < alignment_length - 1:
        start, end = last_end + 1, alignment_length - 1
        if keep(start, end):
            variable_regions.append((start, end))

    return variable_regions

# Function 8: collect_all_variable_regions(...)
# -------------------------------------------------
# Collects variable regions from all *_aligned.fasta files into a CSV.
def collect_all_variable_regions(aligned_dir="aligned_files", block_size=20,
                                 min_length=200, max_length=400,
                                 output_file="all_variable_regions.csv"):
    all_regions = []
    fasta_files = glob(os.path.join(aligned_dir, "*_aligned.fasta"))
    if not fasta_files:
        print("No aligned fasta files found.")
        return

    for file_path in fasta_files:
        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            continue

        sequences = [str(rec.seq).upper() for rec in records]
        min_len = min(len(seq) for seq in sequences)
        sequences = [seq[:min_len] for seq in sequences]

        perfect_indices = get_perfect_match_indices(sequences, window_size=block_size)
        variable_regions = find_variable_regions(perfect_indices, min_len,
                                                 kmer_size=block_size,
                                                 min_length=min_length,
                                                 max_length=max_length)
        for start, end in variable_regions:
            all_regions.append({
                "filename": os.path.basename(file_path),
                "start": start,
                "end": end,
                "length": end - start + 1
            })

    if not all_regions:
        print("No variable regions found.")
        return

    # Save regions to CSV
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["filename", "start", "end", "length"])
        writer.writeheader()
        writer.writerows(all_regions)

    print(f"Saved {len(all_regions)} variable regions to {output_file}")

# Function 9: gc_content(seq)
# -------------------------------------------------
# Calculates GC content (%) of a DNA sequence.
def gc_content(seq):
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0

# Function 10: get_conserved_20mer(seqs, pos, direction)
# -------------------------------------------------
# Finds a conserved 20-mer upstream or downstream of a given position.
# Only returns the sequence if perfectly conserved and gap-free.
def get_conserved_20mer(seqs, pos, direction):
    step = -1 if direction == 'upstream' else 1
    for offset in range(0, 100):  # Search up to 100 bp away
        start = pos + step * offset
        end = start + step * 20
        if step == -1:
            start, end = end, start  # Reverse for slicing
        if start < 0 or end > len(seqs[0]):
            continue
        block = [s[start:end] for s in seqs]
        if all(seq == block[0] and '-' not in seq for seq in block):
            return block[0]
    return None

# Function 11: gc_clamp(seq)
# -------------------------------------------------
# Returns True if at least one of the last two bases is G or C (GC clamp).
def gc_clamp(seq):
    return any(base in 'GC' for base in seq[-2:]) if isinstance(seq, str) else False

# Function 12: wallace_tm(seq)
# -------------------------------------------------
# Calculates melting temperature (Tm) using the Wallace rule.
def wallace_tm(seq):
    seq = seq.upper()
    return 2 * (seq.count('A') + seq.count('T')) + 4 * (seq.count('G') + seq.count('C'))

# Function 13: gc_diff_from_50(gc)
# -------------------------------------------------
# Returns absolute difference from 50% GC content. NaN gets a penalty of 100.
def gc_diff_from_50(gc):
    return abs(gc - 50) if pd.notna(gc) else 100

# Function 14: shannon_entropy(column)
# -------------------------------------------------
# Calculates Shannon entropy of bases in a sequence column.
# High entropy = high diversity of bases.
def shannon_entropy(column):
    counts = Counter(column)
    total = sum(counts.values())
    return -sum((count/total) * math.log2(count/total) for count in counts.values() if count > 0)

# Function 15: extract_annotations(fasta_file)
# -------------------------------------------------
# Extracts metadata from FASTA record descriptions for phylogenetic analysis.
# Assumes descriptions have fields like: >ID|gene_cluster:XYZ|genome_name:ABC|gene_callers_id:123
def extract_annotations(fasta_file):
    rows = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        fields = record.description.split("|")
        entry = {
            "entry_id": fields[0].replace(">", ""),
            "gene_cluster": fields[1].split(":")[1],
            "genome_name": fields[2].split(":")[1],
            "gene_callers_id": fields[3].split(":")[1]
        }
        rows.append(entry)
    return pd.DataFrame(rows)
