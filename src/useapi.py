from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation

from compareseq import compare_sequences, read_sequence
from utils import find_cds

Entrez.email = "yanis.ahnia.etu@univ-lille.fr"
recordGB = Entrez.efetch(db="nucleotide",  id="NM_007389", rettype="gb", retmode="text")
recordFasta = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="fasta", retmode="text")
a=(SeqIO.read(recordGB,'genbank'))
b=(SeqIO.read(recordFasta,'fasta'))

aa=(a.seq)
bb=(b.seq)

print(compare_sequences(aa,bb))

recordGB.close()
recordFasta.close()
print(find_cds(a))

pmid="NM_007389"
handle=Entrez.elink(dbfrom="nucleotide", id=pmid,linkname="pubmed_pubmed")
record = Entrez.read(handle)
print(record[0]["LinkSetDb"][0]["LinkName"])
linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
print (linked)


### Esummary 
handles=Entrez.esummary(db='nucleotide',id='11435')
recordsum=Entrez.read(handles)
handles.close()
print(recordsum[0]['ChrAccVer'])