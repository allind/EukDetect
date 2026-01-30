#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):

	lens = {seq.id: len(str(seq.seq)) for seq in SeqIO.parse(sys.argv[1], 'fasta')}
	taxids = {}
	genomes = {}
	for line in open(sys.argv[2]):
		line = line.strip('\n')
		taxid = line.split('\t')[1]
		if taxid not in taxids:
			taxids[taxid] = []
		taxids[taxid].append(line.split('\t')[0])

		genome = line.split('\t')[2]
		if genome not in genomes:
			genomes[genome] = []
		genomes[genome].append(line.split('\t')[0])

	
	for t in taxids:
		totallen = 0
		for g in taxids[t]:
			totallen += lens[g]
		print(t + '\t' + str(totallen))
	for g in genomes:
		totallen = 0
		for g in genomes[g]:
			totallen += lens[g]
		print(g + '\t' + str(totallen))
if __name__ == "__main__":
  main(sys.argv)
