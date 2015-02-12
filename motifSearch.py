##########
# Usage:
# python motifSearch.py <sequences.fa> <motifs.txt>
##########
import Bio
from Bio import SeqIO
import sys
import MOODS._cmodule
import math

fasta = SeqIO.parse(sys.argv[1], 'fasta')
sequences = {}
for sequence in fasta:
	sequences[str(sequence.id)] = str(sequence.seq)

motifs = {}
motifsfile = open(sys.argv[2], 'r')
for line in motifsfile:
	if line[0] == '>':
		motif_name = line[1:-1]
		motifs[motif_name] = []
	else:
		base = [int(i) for i in line.split()]
		motifs[motif_name].append(base)

motif_list = [motifs[i] for i in motifs]

###Running MOODS
for sequence in sequences:
	print sequence
	for motif in motifs:
		hits = MOODS.search(sequences[sequence], [motifs[motif]], 0.001)
		print 'Motif ' + motif + ': ' + str(len(hits[0]))
	print ''
#END