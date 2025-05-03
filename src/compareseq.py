import argparse
from Bio import SeqIO

def read_sequence(file_path):
    try:
        record = next(SeqIO.parse(file_path, "fasta"))
    except StopIteration:
        record = next(SeqIO.parse(file_path, "genbank"))
    return record.seq

def compare_sequences(seq1, seq2):
    if seq1 == seq2:
        return "The sequences are identical."
    
    min_length = min(len(seq1), len(seq2))
    for i in range(min_length):
        if seq1[i] != seq2[i]:
            return f"The sequences differ at position {i + 1} (1-based index)."
    
    return f"The sequences differ in length: {len(seq1)} vs {len(seq2)}."

def main():
    parser = argparse.ArgumentParser(description="Compare two sequences from FASTA or GenBank files.")
    parser.add_argument("-a", required=True, help="First sequence file (FASTA or GenBank)")
    parser.add_argument("-b", required=True, help="Second sequence file (FASTA or GenBank)")
    
    args = parser.parse_args()
    
    try:
        seq1 = read_sequence(args.a)
        seq2 = read_sequence(args.b)
        result = compare_sequences(seq1, seq2)
        print(result)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
