import argparse

def main():
    # Création du parseur d'arguments
    parser = argparse.ArgumentParser(description="Analyse de promoteurs et recherche de motifs")
    
    # Ajout des options
    parser.add_argument("-t", "--threshold", type=float, required=True,
                       help="Seuil de score")
    parser.add_argument("-l", "--promotor-length", type=int, default=1000,
                       help="Longueur du promoteur (par défaut: 1000)")
    parser.add_argument("-w", "--window-size", type=int, default=40,
                       help="Longueur de la fenêtre glissante (par défaut: 40)")
    
    # Arguments positionnels (liste des mRNA)
    parser.add_argument("mrna_ids", nargs="+", help="Identifiants des mRNA à analyser")
    
    # Analyse des arguments
    args = parser.parse_args()
    
    # Affichage des valeurs
    print("Seuil de score:", args.threshold)
    print("Longueur du promoteur:", args.promotor_length)
    print("Taille de la fenêtre:", args.window_size)
    print("Liste des mRNA:", args.mrna_ids)

if __name__ == "__main__":
    main()