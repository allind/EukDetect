#! /usr/bin/env python
#usage: script.py cluster_file fasta
from Bio import SeqIO
import sys,re
def main(argv):

    #to_replace = {line.strip('\n'): "" for line in open(sys.argv[1])}
    to_replace = {line.split('\t')[0]: "" for line in open(sys.argv[1])}
    to_remove = []
    groups = {}
    seen = {}
    corr = {}
    for line in open(sys.argv[1]):
        genes = line.split('\t')[3].split(',')
        rep = line.split('\t')[1]
        for m in genes:
            if m != rep:
                to_remove.append(m)
        line = line.strip('\n')
        realname = line.split('\t')[1]
        #nothing needed for within genome, just remove the non remove
        name = realname
        if "Same species" in line:
            num = line.split('\t')[2]
            #sp = line.split('\t')[2].split("-")[0]
            sp = re.split('-\d*at\d*-', '-'.join(line.split('\t')[1].split('-')[1:]))[0]
            group = line.split('\t')[1].split('-')[0]
            name = group + '-' + sp + "_SSCollapse_SP" + str(num)
            if name in seen:
                oriname = name
                name = name + '-' + str(seen[name])
                seen[oriname] += 1

            else:
                seen[name] = 1
            to_replace[realname] = name
        elif "Same genus" in line:
            num = line.split('\t')[2]
            genus = line.split('\t')[1].split("_")[0]
            name = genus + "_SPCollapse_SP" + str(num)
            if name in seen:
                oriname = name
                name = name + '-' + str(seen[name])
                seen[oriname] += 1
            else:
                seen[name] = 1
            to_replace[realname] = name
        elif 'Multi genus' in line:
            num = line.split('\t')[2]
            gs = []
            group = line.split('\t')[1].split('-')[0]
            name = group + "_MGCollapse_EASP" + str(num)
            if name in seen:
                oriname = name
                name = name + '-' + str(seen[name])
                seen[oriname] += 1
            else:
                seen[name] = 1

            to_replace[realname] = name
        corr[name] = []
        for g in genes:
            corr[name].append(line)
    correspondence = {}
    #outdest = open(sys.argv[4], 'w')
    for seq in SeqIO.parse(sys.argv[2], 'fasta'):
        if seq.id in to_replace:
            #print(seq.id + '\t' + to_replace[seq.id])
            #corresponence.append(seq.id + '\t' + to_replace[seq.id])
            correspondence[seq.id] = to_replace[seq.id]
            #outdest.write(">" + to_replace[seq.id] + '\n' + str(seq.seq) + '\n')
            print(">" + to_replace[seq.id] + '\n' + str(seq.seq))
        elif seq.id not in to_remove:
            print(">" + seq.id + '\n' + str(seq.seq))
            #outdest.write(">" + seq.id + '\n' + str(seq.seq) + '\n')
    #outdest.close()


    outdest = []
    for line in open(sys.argv[1]):
        line = line.strip('\n')
        ori = line.split('\t')[1].strip('>')
        if ori in correspondence:
            outdest.append(correspondence[ori] + '\t' + line + '\n')
        #outdest.append(line + '\n')
    #outdest.close()
    prefix=sys.argv[2].split('.')[0]
    outfile = open(prefix + "_mmseqs_collapsed_seqnames.txt", 'w')
    for line in outdest:
        outfile.write(line)
    outfile.close()

if __name__ == "__main__":
  main(sys.argv)
