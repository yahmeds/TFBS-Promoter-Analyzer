#!/usr/bin/env python3
import argparse
import json
import numpy as np
from Bio import motifs, Entrez
from Bio.Seq import Seq

from utils import mrna_to_gene, upstream_gene_seq
from pwm import pwm2pssm, scan_all_sequences, best_window, score_window

def json_serializable(obj):
    """
    Convertit des objets non sérialisables en JSON en types compatibles JSON.
    
    Args:
        obj: Objet en entrée à convertir
    
    Returns:
        Une représentation JSON-sérialisable de l'objet
    """
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

def parse_arguments():
     # Création du parseur d'arguments
    parser = argparse.ArgumentParser(description="Analyse de promoteurs et recherche de motifs")
    # Ajout des options
    parser.add_argument("-m", "--pfm", type=str, required=True, help="Fichier contenant la matrice")
    parser.add_argument("mrna_ids", nargs="+", help="Identifiants des mRNA à analyser")
    parser.add_argument("-t", "--threshold", type=float, required=True, help="Seuil de score")
    parser.add_argument("-l", "--promotor-length", type=int, default=1000,
                       help="Longueur du promoteur (par défaut: 1000)")
    parser.add_argument("-w", "--window-size", type=int, default=40,
                       help="Longueur de la fenêtre glissante (par défaut: 40)")
    parser.add_argument("-s", "--window-threshold", type=float,
                       help="Seuil de score de fenetre")
    parser.add_argument("-p", "--pseudocount", type=float, default=0.1,
                       help="Valeur du pseudo-poids")
    parser.add_argument("-o", "--output", help="Nom du fichier JSON des résultats") 
    
    args = parser.parse_args()
    
    return args

def get_promoter_sequences(mrna_ids, promoter_length):
    """
        Récupère les séquences promoteurs pour des identifiants d'ARNm donnés.
        
        Args:
            mrna_ids (list): Liste d'identifiants d'ARNm
            promoter_length (int): Longueur de la séquence promoteur à récupérer
        
        Returns:
            list: Liste des séquences promoteurs
    """
    Entrez.email = "m-hand.ait-cherif.etu@univ-lille.fr" 
    promoter_sequences = []
    
    for mrna_id in mrna_ids:
        try:
            gene_id = mrna_to_gene(mrna_id)
            seq = upstream_gene_seq(gene_id, promoter_length)
            seq.id = mrna_id
            promoter_sequences.append(seq)
        except Exception as e:
            print(f"Error processing {mrna_id}: {e}")
    
    return promoter_sequences

def process_tfbs(args, promoter_sequences):
    """
        Fonction principale de traitement pour l'identification de sites de liaison 
        putatifs aux facteurs de transcription (TFBS).
        
        Args:
            args (argparse.Namespace): Arguments passés en ligne de commande (par argparse)
            promoter_sequences (list): Liste des séquences promoteurs
        
        Returns:
            dict: Résultats des TFBS traités dans un format structuré
    """
    with open(args.pfm) as handle:
        matrix = motifs.read(handle, "jaspar")
        pssm = pwm2pssm(matrix.counts, args.pseudocount)
    
    scan_results = scan_all_sequences(pssm, promoter_sequences, args.threshold)
    # création dictinnaire results 
    results_dict = {
        "parameters": {
            "jaspar_file": args.pfm,
            "threshold": json_serializable(args.threshold),
            "window_size": json_serializable(args.window_size),
            "window_threshold": json_serializable(args.window_threshold),
            "promoter_length": json_serializable(args.promotor_length),
            "pseudocount": json_serializable(args.pseudocount),
            "mrna_ids": args.mrna_ids
        },
        "matches": [],
        "windows": []
    }    
    all_matches = []
    for i, seq_scan_results in enumerate(scan_results):
        sorted_results = sorted(seq_scan_results, key=lambda x: x[0])
        if not sorted_results:
            continue
        for pos, score in sorted_results:
            if score >= args.threshold:
                match = {
                    "position": json_serializable(pos),
                    "score": json_serializable(score),
                    "mrna_id": promoter_sequences[i].id,
                    "motif": matrix.name
                }
                all_matches.append(match)
                results_dict["matches"].append(match)
    
    all_matches.sort(key=lambda x: abs(x['position']))    
    
    def find_windows(matches, window_size):
        """ 
            Trouve des fenêtres de taille constante en fonction de la taille de fenêtre spécifiée, en veillant à ce que toutes les fenêtres soient de la même longueur et positionnées de manière optimale.

            Args : 
                matches (list) : Liste de dictionnaires de correspondances 
                window_size (int) : Taille de fenêtre souhaitée

            Returns : 
                list : Liste de fenêtres, chacune contenant des motifs correspondants 
        """
        windows = []
        if not matches:
            return windows
        sorted_matches = sorted(matches, key=lambda x: abs(x['position']))
        while sorted_matches:
            seed_match = sorted_matches.pop(0)
            seed_pos = seed_match['position']
            
            window_start = seed_pos - window_size // 2
            window_end = seed_pos + window_size // 2
        
            window_matches = [seed_match]
            
            for match in sorted_matches[:]:
                if window_start <= match['position'] <= window_end:
                    window_matches.append(match)
                    sorted_matches.remove(match)
            
            windows.append(window_matches)
        
        return windows
    
    windows = find_windows(all_matches, args.window_size)
    
    for window_num, window in enumerate(windows, 1):
        window_start = min(m['position'] for m in window) - (args.window_size // 2)
        window_end = max(m['position'] for m in window) + (args.window_size // 2)

        for match in window:
            match['window_id'] = window_num
        
        start = window_start
        end = window_end

        updated_score = score_window(scan_results, start, end)
        
        window_details = {
            "number": window_num,
            "motif": window[0]['motif'],
            "range": f"[{window_start}:{window_end}]",
            "score": updated_score,
            "match_count": len(window),
            "matches": window
        }
        
        results_dict["windows"].append(window_details)
    
    if args.output:
        output_file = args.output if args.output.endswith('.json') else args.output + '.json'
        with open(output_file, 'w') as outfile:
            json.dump(results_dict, outfile, indent=2, default=json_serializable)
    
    return results_dict

def format_results(results, window_size=30):
    """
        Formate les résultats de TFBS en une table structurée.
        
        Args:
            results (dict): Résultats traités des TFBS
    """
    results['matches'].sort(key=lambda x: abs(x['position']))
    
    for window in results['windows']:
        print(f"{window['number']} {window['motif']} {window['range']} {window['score']:.2f} {window['match_count']}")
    print() 
    for window in results['windows']:
        print()
        sorted_matches = sorted(window['matches'], key=lambda m: m['position'])
        for match in sorted_matches:
            print(f"{window['number']} {match['mrna_id']} {match['motif']} {match['position']} {match['score']:.2f}")

def main():
    args = parse_arguments()
    promoter_sequences = get_promoter_sequences(args.mrna_ids, args.promotor_length)
    results = process_tfbs(args, promoter_sequences)
    format_results(results, args.window_size)
    if not args.output:
        print(json.dumps(results, indent=2, default=json_serializable))

if __name__ == "__main__":
    main()