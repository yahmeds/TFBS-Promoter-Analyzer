from Bio import motifs
from Bio.Seq import Seq
from Bio import Entrez
from utils import mrna_to_gene, upstream_gene_seq




def scan_all_sequences(pssm, sequences, threshold):
    """
    Recherche les occurrences d'une PSSM dans une liste de séquences.
    Retourne un dictionnaire où chaque clé est une séquence et la valeur
    est une liste de tuples (position, score) pour les occurrences.

    Args:
        pssm (PositionSpecificScoringMatrix): La PSSM utilisée pour la recherche.
        sequences (list of Bio.Seq): Liste des séquences à analyser.
        threshold (float): Seuil de score minimal pour considérer une occurrence.

    Returns:
        list: Liste des résultats de scan, chaque élément étant une liste de tuples (position, score).
              Chaque sous-liste correspond à une séquence de la liste `sequences`.
    """
    scan_results = []
    for seq in sequences:
        # Pour chaque séquence nous recherchons les occ de la PSSM et appliquons scan_sequence
        scan_results.append(scan_sequence(pssm, seq, threshold))
    return scan_results


def pwm2pssm(freq_matrix, pseudocount=0.1):

    """
    Convertit une matrice de fréquence en PSSM avec pseudo-comptes

    Args:
        freq_matrix (FrequencyPositionMatrix): Matrice de fréquence JASPAR
        pseudocount (float): Valeur des pseudo-poids, par défaut 0.1
        
    Returns:
        PositionSpecificScoringMatrix: La PSSM correspondante
    """
    #Normalisation de la matrice de fréquence avec les pseudo count
    pwm = freq_matrix.normalize(pseudocount)
    #Transformation en PSSM via log_odds
    pssm = pwm.log_odds()
    
    return pssm

    

def scan_sequence(pssm, sequence, threshold):
    """
    Recherche les occurrences d'une PSSM dans une séquence avec un seuil donné
    Retourne une liste de tuples (position, score) pour le brin positif uniquement
    
    Args:
        pssm (PositionSpecificScoringMatrix): La PSSM à rechercher
        sequence (Bio.Seq): La séquence à analyser
        threshold (float): Le seuil de score minimal
        
    Returns:
        list: Liste de tuples (position, score) des occurrences   
             Chaque tuple contient une position et un score d'occurrence. 
    """
    # Recherche uniquement sur le brin positif (strand='+')
    occurrences = []
    length = len(sequence)
    for position, score in pssm.search(sequence, threshold=threshold,both= True):
        # La position est relative au début de la séquence promotrice
        # transformation en pos en amont du gene 
        aposition = int(position - length) 
        
        if aposition >= -length:
            # Si la position est valide dans le cadre de la séquence promotrice nous 
         #   ajoutons l'occurrence sous forme (position, score)
           
            occurrences.append((aposition, float(score)))  
    return occurrences

def score_window(scan_results, start, end):
    """
    Calcule le score d'une fenêtre à partir des résultats de scan_all_sequences
    
    Args:
        scan_results (list): Résultat de scan_all_sequences, liste contenant des sous-listes de tuples (position, score).
        start (int): Position de début de la fenêtre
        end (int): Position de fin de la fenêtre
        
    Returns:
        float: Score cumulé pondéré des occurrences dans la fenêtre.
    """
    total_score = 0.0
    
    for occurrences in scan_results:
        for pos, score in occurrences:
            if start <= pos <= end:
                # Poids basé sur l'intensité du score
                poids = abs(score)  
                 # Ajout du score pondéré à la somme totale
                total_score += score * poids
    
    return total_score



def best_window(scan_results, window_size):
    """
    Trouve la fenêtre de meilleure score
     
    Args:
        scan_results (list): Résultat de scan_all_sequences

        window_size (int): taille de la fenêtre        
    Returns:
        float: meilleure position de debut pour la fenêtre
        float: meilleure position de fin  pour fenêtre
        float: meilleure score  des occurrences dans la fenêtre
    """
    
    best_start, best_end = None, None
  #meilleur score initalisé a une petite valeur    
    best_score = float('-inf')

    all_occurrences = []
    
    for occurrences in scan_results:
       #ici on extrait uniquement les positions  
        all_occurrences.extend([pos for pos, _ in occurrences])

    if not all_occurrences:
        return None, None, None
   # on trie les positions afin de parcourir les fenetres dans l'ordre en supprimant les doublons s'ils existent
    all_occurrences = sorted(set(all_occurrences))


   # création de la fenetre 
    for occ in all_occurrences:
       #on place pos au milieu de notre fenetre  
        current_occ = int(occ)  # Convert to Python int
        start = current_occ - window_size // 2
        end = current_occ + window_size // 2
        current_score = score_window(scan_results, start, end)      

        print(f"Fenêtre autour de l'occurrence à {occ}: Début = {start}, Fin = {end}, Score = {current_score}")
        if current_score > best_score:
            best_score = current_score
            best_start, best_end = start, end


    
    return best_start, best_end, best_score

if __name__ == "__main__":
    Entrez.email = "yanis.ahnia.etu@univ-lille.fr"
    length = 1000
    liste =["NM_000451", "NM_007389" ]
    file = "./data/MA0083.3.jaspar"
    threshold = -20
    l= []
    for id in liste : 
        gene_id = mrna_to_gene(id)
        seq = upstream_gene_seq(gene_id,length)
        seq.id = id
        l.append(seq)

    with open(file) as handle:
        matrix = motifs.read(handle, "jaspar")
        pssm = pwm2pssm(matrix.counts,0.1)

        results = scan_all_sequences(pssm,l,threshold)
        print("SCAN DES SEQUENCES:", results)
        # Calculer le score pour une fenêtre spécifiée
        score_result = score_window(results, -200, -100)  
        print(f"Score de la fenêtre: {score_result}")
        
        # Recherche de la meilleure fenêtre
        best_window_result = best_window(results, window_size=40)  
        print("Meilleure fenêtre trouvée:", best_window_result)
 