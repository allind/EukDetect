#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):
	strings = []
	removed = []
	#fdest = open(sys.argv[1].split('.')[0] + "_rmmask.fasta", 'w')
	for seq in SeqIO.parse(sys.argv[1], 'fasta'):
		length = len(str(seq.seq))
		ns = str(seq.seq).count("N")
		percent = round((ns / length) * 100, 2)
		if percent < 10:
			print(">" + seq.id + '\n' + str(seq.seq) + '\n')
		else:
			removed.append(seq.id)
	
	r = open("removed_for_repetitive_elements.txt", 'w')
	for f in removed:
		r.write(f + '\n')
	
if __name__ == "__main__":
  main(sys.argv)
