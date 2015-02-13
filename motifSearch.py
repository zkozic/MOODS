##########
# Python interface for MOODS package.
##########
import Bio
from Bio import SeqIO
import sys
import MOODS._cmodule
import math
import argparse

##### Parsing command line arguments #####
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sequences', action = 'store', dest = 'sequences', help = 'File with sequences in FASTA format')
parser.add_argument('-m', '--motifs', action = 'store', dest = 'motifs', help = 'File with motifs in JASPAR count format')
parser.add_argument('-t', '--threshold', type = float, action = 'store', dest = 'threshold', default = 0.001, help = 'Threshold value for matrix scanning. Default = 0.001')
parser.add_argument('-o', '--output', action = 'store', dest = 'output', default = './moods_out.txt', help = 'Result output file. Outputs the results in moods_out.txt in current working directory')
parser.add_argument('-b', '--background', action = 'store', dest = 'background', default = None, help = 'Background distribution as an array of four doubles, corresponding to the frequencies of A, C, G and T, respectively. By default the background is estimated from the sequence.')
parser.add_argument('-l', '--log_base', action = 'store', dest = 'logbase', default = None, help = 'Base for logarithms used in log-odds computations. Relevant if using convert_log_odds=True and threshold_from_p=False. Defaults to natural logarithm if None is given.')
parser.add_argument('-p', '--pseudocount', action = 'store', type = int, default = 1, dest = 'pcount', help = 'Pseudocount used in log-odds conversion and added to sequence symbol counts when estimating the background from sequence. Default 1')
parser.add_argument('-B', '--both_strands', action = 'store_true', dest = 'strands', default = False, help = 'Scans against reverse complement sequence in addition to the input sequence. Hits on reverse complement are reported at position [position - sequence_length], which is always negative. The actual hit site for any hit is always seq[pos, pos + matrix_length]. Default False.')
args = parser.parse_args()
##########################################

fasta = SeqIO.parse(args.sequences, 'fasta')
sequences = {}
for sequence in fasta:
	sequences[str(sequence.id)] = str(sequence.seq)

motifs = {}
motifsfile = open(args.motifs, 'r')
for line in motifsfile:
	if line[0] == '>':
		motif_name = line[1:-1]
		motifs[motif_name] = []
	else:
		base = [int(i) for i in line.split()]
		motifs[motif_name].append(base)

motif_list = [motifs[i] for i in motifs]

###Running MOODS
output = open(args.output, 'w')
for sequence in sequences:
	output.write(sequence + '\n')
	for motif in motifs:
		hits = MOODS.search(sequences[sequence], [motifs[motif]], thresholds = args.threshold, bg = args.background, both_strands = args.strands, log_base = args.logbase, pseudocount = args.pcount)
		output.write('Motif ' + motif + ': ' + str(len(hits[0])) + '\n')
	output.write('\n')
#END