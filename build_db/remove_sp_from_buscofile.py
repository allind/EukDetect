#! /usr/bin/env python
#usage: script.py [busco fasta] [list of species to remove]

from Bio import SeqIO
import sys,re
def main(argv):
    
    remove = [line.strip('\n') for line in open(sys.argv[2])]
    for seq in SeqIO.parse(sys.argv[1], 'fasta'):
        sp = re.split('-\d*at\d*-', '-'.join(seq.id.split('-')[1:]))[0]
        if sp not in remove:
            print(">" + seq.id)
            print(str(seq.seq))

if __name__ == "__main__":
  main(sys.argv)
