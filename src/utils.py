from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio import Entrez
from Bio.Seq import Seq
import os 


def find_cds(seq_record):
    """
    Trouve les positions des CDS dans un SeqRecord
    
    Args:
        record (SeqRecord): Le SeqRecord à analyser
        
    Returns:
        list: Liste de tuples (start, end) des CDS
    """
    cds_positions = []
    
    for feature in seq_record.features:
        if feature.type == "CDS": 
            start = int(feature.location.start)
            end = int(feature.location.end)
            cds_positions.append((start, end))
    
    return cds_positions

def mrna_to_gene(mrna_accession):
    """
    Trouve l'identifiant du gène correspondant à un ARNm
    
    Args:
        mrna_accession (str): Numéro d'accession de l'ARNm
        
    Returns:
        str: Identifiant du gène
    """
    Entrez.email = "yanis.ahnia.etu@univ-lille.fr"
    handle = Entrez.elink(dbfrom="nucleotide", id=mrna_accession, db="gene")
    result = Entrez.read(handle)
    handle.close()
    gene_id = result[0]["LinkSetDb"][0]["Link"][0]["Id"]
    return gene_id

def upstream_gene_seq(gene_id, length):
    """
    Récupère la séquence en amont d'un gène
    
    Args:
        gene_id (str): Identifiant du gène
        length (int): Longueur de la séquence à récupérer
        
    Returns:
        Bio.Seq: Séquence en amont du gène
    """
    # Récupérer les informations sur le gène
    Entrez.email = "yanis.ahnia.etu@univ-lille.fr"
    handle = Entrez.esummary(db="gene", id=gene_id)
    result = Entrez.read(handle)
    handle.close()
    
    genomic_info = result["DocumentSummarySet"]["DocumentSummary"][0].get("GenomicInfo", [])
    if genomic_info:
        chrom_accession = genomic_info[0].get("ChrAccVer")
        start = int(genomic_info[0].get("ChrStart", 0))
        end = int(genomic_info[0]["ChrStop"])           # Position de fin du gène
        strand = int(genomic_info[0].get("Strand", 1))  # Ajout pour récupérer le brin

    if start < end :
        endPos = start - length +1
        handle1 = Entrez.efetch(db="nucleotide", id = chrom_accession, rettype = "gb", retmode ="texte", seq_start=str(start), seq_stop=str(endPos), strand ="1")

    else: 
        endPos = start + length +1
        handle1 = Entrez.efetch(db="nucleotide", id = chrom_accession, rettype = "gb", retmode ="texte", seq_start=str(start+2), seq_stop=str(endPos), strand ="2")
    

    record1 = SeqIO.read(handle1, 'gb')
    sequence = record1.seq

    return sequence

#cas en amont du promotteur

def upstream_gene_seq2(gene_id, length):
    """
    Récupère la séquence en amont d'un gène
    
    Args:
        gene_id (str): Identifiant du gène
        length (int): Longueur de la séquence à récupérer
        
    Returns:
        Bio.Seq: Séquence en amont du gène
    """
    # Récupérer les informations sur le gène
    Entrez.email = "yanis.ahnia.etu@univ-lille.fr"
    handle = Entrez.esummary(db="gene", id=gene_id)
    result = Entrez.read(handle)
    handle.close()
    
    genomic_info = result["DocumentSummarySet"]["DocumentSummary"][0].get("GenomicInfo", [])
    if genomic_info:
        chrom_accession = genomic_info[0].get("ChrAccVer")
        start = int(genomic_info[0].get("ChrStart", 0))
        strand = int(genomic_info[0].get("Strand", 1))  # Ajout pour récupérer le brin

    if strand == 1:
        seq_start = max(0, start - length + 1)  # Prendre en amont pour le brin positif
        seq_end = start

    else:
        seq_start = start
        seq_end = min(start + length, start + length) 

        
    handle = Entrez.efetch(db="nucleotide", id=chrom_accession, rettype="fasta", retmode="text", seq_start=seq_start, seq_stop=seq_end)
    record = SeqIO.read(handle, "fasta")
    handle.close()


    return record.seq if strand == 1 else record.seq.reverse_complement() 

def download_promotors(id_list, length, output_dir="."):
    """
    Télécharge les séquences promotrices d'une liste d'ARNm
    
    Args:
        mrna_list (list): Liste d'identifiants d'ARNm
        promotor_length (int): Longueur des séquences promotrices
        output_dir (str): Répertoire de sortie, par défaut '.'
    """
    for i in id_list:
        try:
            gene_id = mrna_to_gene(i)
            promoter_seq = upstream_gene_seq(gene_id, length)
            file_path=os.path.join(output_dir, f"{i}_{length}.fa")
            with open(file_path, "w") as f:
                f.write(f">{i}_promoter_{length}\n{str(promoter_seq)}\n")
            print(f"Fichier enregistré : {file_path}")
        except Exception as e:
            print(f"Erreur lors du téléchargement pour {i}: {e}")




def main():
    record1 = SeqIO.read("data/sequence.gb", "genbank")
    record2 = SeqIO.read("data/sequence_NM_000451.4.gb", "genbank")
    cds1 = find_cds(record1)
    cds2 = find_cds(record2)
    print(f"CDS dans sequence.gb : {cds1}")
    print(f"CDS dans sequence_NM_000451.4.gb : {cds2}")
    gene_id = mrna_to_gene("NM_007389")
    print(f"Identifiant du gène correspondant à NM_007389 : {gene_id}")
    upstream_seq = upstream_gene_seq("11435", 250)  # Exemple avec 1000 pb
    print(f"Séquence en amont du gène {gene_id} :\n{upstream_seq}")

    print("----------------------------------------------------------------------------")  
    download_promotors(["NM_007389", "NM_079420", "NM_001267550", "NM_002470", "NM_003279", "NM_005159", "NM_003281", "NM_002469", "NM_004997", "NM_004320", "NM_001100", "NM_006757"], 1000, "data/promoters")

if __name__ == "__main__":
    main()
    
