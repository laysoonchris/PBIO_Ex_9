import random
import csv
import re
from typing import Optional

try:
    import matplotlib.pyplot
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

def validate_positive_integer(prompt: str, min_value: int=1, max_value: int=100_000) -> int:
    while True:
        raw_input = input(prompt)
        try:
            value = int(raw_input)
            if min_value <= value <= max_value:
                return value
            print(f"Please enter an integer between {min_value} and {max_value}.")
        except ValueError:
            print(f"Please enter an integer between {min_value} and {max_value}.")

def validate_seq_id(prompt: str) -> str:
    while True:
        seq_id = input(prompt).strip()
        if seq_id and not re.search(r"\s", seq_id):
            return seq_id
        print("Please enter a valid sequence ID (without whitespace).")

def generate_random_sequence(length: int, weights: Optional[dict] = None) -> str:
    nucleotides = ['A', 'C', 'G', 'T']
    if weights:
        w = [weights["A"], weights["C"], weights["G"], weights["T"]]
        return "".join(random.choices(nucleotides, weights=w, k=length))
    
    return "".join(random.choices(nucleotides, k=length))

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    header = f">{seq_id} {description}" if description else f">{seq_id}"
    lines = [sequence[i: i+line_width] for i in range(0, len(sequence), line_width)]

    return header + "\n" + "\n".join(lines)

def calculate_stats(sequence: str) -> dict:
    bio_seq = "".join(c for c in sequence if c in "ACGT")
    total_length = len(bio_seq)

    if total_length == 0:
        return {
            "A": 0.0,
            "C": 0.0,
            "G": 0.0,
            "T": 0.0,
            "GC_content": 0.0
        }
    
    counts = {nucleotide: bio_seq.count(nucleotide) for nucleotide in "ACGT"}
    stats = {nucleotide: round(counts[nucleotide] / total_length * 100, 2) for nucleotide in "ACGT"}
    stats["GC"] = round((counts['G'] + counts['C']) / total_length * 100, 2)
    return stats

def print_stats(stats: dict, length: int) -> None:
    print(f"\nStatistics for sequence (length: {length}):")
    for nucleotide in ["A", "C", "G", "T"]:
        print(f"{nucleotide}: {stats[nucleotide]:.2f}%")
    print(f"GC Content: {stats['GC']:.2f}%")

def insert_name(sequence: str, name: str) -> str:
    position = random.randint(0, len(sequence))
    return sequence[:position] + name.lower() + sequence[position:]

# feature 1: batch mode



# feature 2: configurable nucleotide distribution



# feature 5: transcribtion in silico

def transcribe(dna_sequence: str) -> str:
    return dna_sequence.replace("T", "U")


def sliding_window_gc_content(sequence: str, window_size: int) -> list:
    bio_seq = "".join(c for c in sequence if c in "ACGT")
    gc_content_values = []

    for i in range(len(bio_seq) - window_size + 1):
        window = bio_seq[i: i+window_size]
        gc = (window.count("G") + window.count("C")) / window_size * 100
        gc_content_values.append({"start position": i + 1, "GC-content": round(gc, 2)})
    return gc_content_values

def save_gc_content_to_csv(gc_content_values: list, filename: str) -> None:
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=["start_position", "gc_content"])
        writer.writeheader()
        writer.writerows(gc_content_values)
    print(f"GC results from sliding window saved to: {filename}")

# feature 8: GC-content plot


def main():
    print("Welcome to the DNA Sequence Generator!")