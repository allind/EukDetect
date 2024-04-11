#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):
	corr = {line.split('\t')[1].strip('\n'): line.split('\t')[0].strip('\n') for line in open(sys.argv[1])}
	for seq in SeqIO.parse(sys.argv[2], 'fasta'):
		new = corr[seq.id]

		print(">" + new)
		print(str(seq.seq))
	
if __name__ == "__main__":
  main(sys.argv)
