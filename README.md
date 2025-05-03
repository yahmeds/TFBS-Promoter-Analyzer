# Bioinformatics Project

## Final Submission: 28/03/2025

### Project Description

This project is a Python application designed to analyze gene promoter sequences to identify potential **Transcription Factor Binding Sites (TFBS)**.

**Project Objectives:**
- Automate the retrieval of promoter sequences from mRNA identifiers using biological databases (NCBI).
- Identify TFBS using Position Weight Matrices (PWM) from the JASPAR database.
- Detect promoter regions with high gene regulation potential by analyzing associated scores.

---

### Project Structure

bioinfo/
├── data/
│   └── MA0083.3.jaspar   
├── src/
│   ├── utils.py          
│   ├── pwm.py           
│   └── putative_TFBS.py  
│   └── *.py(Pour les tps precedents) 
├── tests/
│   ├── test_utils.py            
│   ├── test_pwm.py              
│   └── test_putative_tfbs.py    
├── Makefile                      
└── README.md 
└── viewPutativeTFBS.html 


---

### Key Functional Components

#### 1. `utils.py`

This module provides helper functions for:

- Converting mRNA identifiers to NCBI gene IDs.
- Automatically retrieving promoter sequences via the NCBI API.
- Downloading promoter sequences in FASTA format.

**Key Functions:**
- `mrna_to_gene(mrna_accession)`
- `upstream_gene_seq(gene_id, length)`
- `download_promotors(id_list, length, output_dir)`

#### 2. `pwm.py`

Handles PWM data (e.g., from JASPAR) and converts them into score matrices (PSSM) to scan promoter sequences.

**Main Functions:**
- `pwm2pssm(freq_matrix)` – Convert PWM to PSSM
- `scan_sequence(pssm, sequence, threshold)` – Scan a single sequence
- `scan_all_sequences(pssm, sequences, threshold)` – Scan multiple sequences
- Evaluation and identification of high-scoring promoter regions

#### 3. `putative_TFBS.py`

Main module that performs the full analysis:

- Automatically extracts promoter sequences from given mRNA identifiers.
- Scans for potential TFBS using a specified matrix.
- Outputs structured results in JSON format.
- Ignores windows with identical scores.

**Key Functions:**
- `get_promoter_sequences(mrna_ids, promoter_length)`
- `process_tfbs(args, promoter_sequences)`
- `format_results(results)`

---

### Example Execution

```bash
python3 src/putativeTFBS.py -m data/MA0083.3.jaspar -t -20 -s 10 NM_000451 NM_007389 -o output


### Exemple d'éxecution:
```bash 
python3 src/putativeTFBS.py -m data/MA0083.3.jaspar -t -20 -s 10 NM_000451 NM_007389 -o output
``` 
### Run tests

```bash
  make test
```
### Visualizing the Generated JSON
When hovering the mouse over an occurrence, an information tooltip appears.
![alt](Affichage_Curseur.png)

## JSON File Structure

### 1. Paramètres (`parameters`)
```json
{
  "parameters": {
    "jaspar_file": "data/MA0083.3.jaspar",
    "threshold": -20.0,
    "window_size": 30,
    "promoter_length": 1000,
    "mrna_ids": [
      "NM_000451", 
      "NM_007389"
    ]
  }
}
```

### 2. Correspondances (`matches`)
Chaque correspondance représente un site de liaison individuel :
```json
{
  "position": -922,
  "score": -15.604,
  "mrna_id": "NM_000451",
  "motif": "SRF",
  "window_id": 17
}
```

### 3. Fenêtres (`windows`)
Regroupement des correspondances par fenêtres :
```json
{
  "number": 17,
  "motif": "SRF",
  "range": "[907:937]",
  "score": -15.604,
  "match_count": 1,
  "matches": [...]
}
```

## 🔬 Visualization Process

The `drawTFBS()` function creates a graphical view based on the generated JSON data.

### 1. Initialization
- Loads the list of genes from the `parameters` section
- Defines the length of the promoter region

### 2. Drawing TFBS Sites
- Iterates through the `matches`
- Positions each TFBS site on a horizontal axis
- Uses score-based color coding (darker blue = higher relevance)

### 3. Drawing Windows
- Iterates through the `windows`
- Locates the first match in each window
- Draws a red rectangle around the window region
- Labels each window with its number

---

## 💻 Usage Instructions

1. Generate your JSON output using the analysis script.
2. Load the file into the provided HTML interface.
3. The visualization will adapt dynamically to the content of the JSON.

---

## 💡 Tips

- **Darker blue** = higher TFBS match score
- **Red rectangles** = important windows with clustered sites
- **Hover over** a site to see detailed info (position, score, motif, etc.)

---

## 🧠 Technical Notes

- TFBS positions are negative, indicating upstream regions relative to the gene start site
- Scores reflect similarity to the input motif
- The visualization supports **multiple genes simultaneously**

---

## 📘 Answers

**2.** Using `record.seq`, the sequences are identical.  
**3.** Use the `reverse_complement()` function  
**4.** Use `.id` to retrieve sequence identifiers  
**5.** Access positions from features using:

```python
start = feature.location.start
end = feature.location.end


### API du NCBI with Biopython
#### efetch
1. The sequences retrieved are indeed identical.
#### elink
2.  Convert mRNA to Gene ID using elink:
```
pmid="NM_007389"
handle=Entrez.elink(dbfrom="nucleotide", id=pmid,linkname="pubmed_pubmed")
record = Entrez.read(handle)
linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
```
