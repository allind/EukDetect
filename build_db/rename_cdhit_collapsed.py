#! /usr/bin/env python
#usage: script.py cluster_file fasta
from Bio import SeqIO
import sys,re
def main(argv):

    #to_replace = {line.strip('\n'): "" for line in open(sys.argv[1])}
    to_replace = {line.split('\t')[0]: "" for line in open(sys.argv[1])}
    to_remove = []
    seen = {}
    for line in open(sys.argv[1]):
        rm = line.split('\t')[3:-1]
        for m in rm:
            to_remove.append(m)
        line = line.strip('\n')
        realname = line.split('\t')[0]
        if "same_species" in line:
            num = len(line.split('\t')[2:])
            #sp = line.split('\t')[2].split("-")[0]
            sp = re.split('-\d*at\d*-', '-'.join(line.split('\t')[2].split('-')[1:]))[0]
            group = line.split('\t')[2].split('-')[0]
            name = group + '-' + sp + "_SSCollapse_SP" + str(num)
            if name in seen:
                oriname = name
                name = name + '-' + str(seen[name])
                seen[oriname] += 1

            else:
                seen[name] = 1
            to_replace[realname] = name
        elif "same_genus" in line:
            num = len(line.split('\t')[2:])
            genus = line.split('\t')[2].split("_")[0]
            name = genus + "_SPCollapse_SP" + str(num)
            if name in seen:
                oriname = name
                name = name + '-' + str(seen[name])
                seen[oriname] += 1
            else:
                seen[name] = 1
            to_replace[realname] = name
        else:
            num = len(line.split('\t')[2:])
            gs = []
            for e in line.split('\t')[2:]:
                g = e.split('_')
                if g not in gs:
                    gs.append(g)
                name = "_".join(g) + "_MGCollapse_EASP" + str(num)
                if name in seen:
                    oriname = name
                    name = name + '-' + str(seen[name])
                    seen[oriname] += 1
                else:
                    seen[name] = 1

                to_replace[realname] = name
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

    #outdest = open("busco_collapsed_sequence_id_changes_v4.txt", 'w')
    #foar e in correspondence:
#        outdest.write(e + '\n')
#    outdest.close()

    outdest = []
    for line in open(sys.argv[1]):
        line = line.strip('\n')
        ori = line.split('\t')[0].strip('>')
        outdest.append(correspondence[ori] + '\t' + line + '\n')
        #outdest.append(line + '\n')
    #outdest.close()
    prefix=sys.argv[2].split('.')[0]
    outfile = open(prefix + "_cdhit_collapsed_seqnames.txt", 'w')
    for line in outdest:
        outfile.write(line)
    outfile.close()

if __name__ == "__main__":
  main(sys.argv)
