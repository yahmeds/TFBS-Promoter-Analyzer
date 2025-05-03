#!/usr/bin/env python3
import argparse
from Bio import motifs
from Bio.Seq import Seq
import sys
from utils import mrna_to_gene, upstream_gene_seq
from pwm import pwm2pssm,scan_sequence,score_window,best_window,scan_all_sequences

def main():
    parser = argparse.ArgumentParser(description="Scan promoter sequence for PWM matches")
    parser.add_argument("jaspar_file", help="JASPAR matrix file")
    parser.add_argument("mrna_id", help="mRNA GenBank ID")
    parser.add_argument("upstream_length", type=int, help="Length of upstream sequence")
    parser.add_argument("threshold", type=float, help="Score threshold")
    args = parser.parse_args()
    
    try:
        # Lecture des matrices JASPAR
        with open(args.jaspar_file) as handle:
            jaspar_motifs = motifs.parse(handle, "jaspar")
        
        # Récupération de la séquence promotrice
        gene_id = mrna_to_gene(args.mrna_id)
        upstream_seq = upstream_gene_seq(gene_id, args.upstream_length)
        scan_results = {}
        # Pour chaque matrice dans le fichier JASPAR
        for motif in jaspar_motifs:
            # Conversion en PSSM
            pssm = pwm2pssm(motif.counts, pseudocount=0.1)
            
            # Recherche des occurrences
            occurrences = scan_sequence(pssm, upstream_seq, args.threshold)
            scan_results[motif.name] = occurrences  
                        
            # Affichage des résultats
            for pos, score in occurrences:
                print(f"{motif.name} {pos} {score}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
