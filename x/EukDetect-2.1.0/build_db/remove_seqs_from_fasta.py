#! /usr/bin/env python
#usage: script.py [to_remove] [fasta]
from Bio import SeqIO
import sys
def main(argv):

	remove = [line.strip('\n') for line in open(sys.argv[1])]

	for seq in SeqIO.parse(sys.argv[2], 'fasta'):
		if seq.id not in remove:
			print(">" + seq.id)
			print(str(seq.seq))
if __name__ == "__main__":
  main(sys.argv)
