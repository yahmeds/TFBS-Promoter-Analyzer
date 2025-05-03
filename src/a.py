from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from Bio.SeqFeature import SeqFeature 

from compareseq import compare_sequences, read_sequence

record1 = SeqIO.read("data/sequence.gb","genbank")
record2 = SeqIO.read("data/sequence.fasta","fasta")
a=(record1.seq)
b=(record2.seq)

record3 = SeqIO.read("data/sequence_NM_000451.4.gb","genbank")
record4 = SeqIO.read("data/sequence_NM_000451.4.fasta","fasta")
c=(record1.seq)
d=(record2.seq)

#print(compare_sequences(a,b))
#print(compare_sequences(c,d))

#print(record1.reverse_complement())
#print(record2.id)
#print(record1.id)
print(SeqFeature(record1.features))
