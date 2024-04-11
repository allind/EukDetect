#! /usr/bin/env python
from Bio import SeqIO
import sys
def main(argv):

    counter = 0
    length = 9
    corr = {}
    for seq in SeqIO.parse(sys.argv[1], 'fasta'):
        counter += 1
        name = (length-len(str(counter)))*'0' + str(counter)
        print(">" + name)
        print(str(seq.seq))
        corr[seq.id] = name
    prefix=sys.argv[1].split('.')[0]
    out = open(prefix + "_busco_seqid_sequential_correspondence.txt", 'w')
    for c in corr:
        out.write(c + '\t' + corr[c] + "\n")
    out.close()
if __name__ == "__main__":
  main(sys.argv)
