#! /usr/bin/env python
import sys
def main(argv):
    

    for line in open(sys.argv[1]):
        line = line.strip('\n')
        dest = open(line.split('\t')[-1].replace(".fasta", "_table.bed"), 'w')
        prefix = line.split('\t')[1]
        spname = line.split('\t')[0]
        table = line.split('\t')[3]
        genome = line.split('\t')[2]
        seen = {}
        for line in open(table):
            line = line.strip('\n')
            if "#" not in line and ("Complete" in line or "Duplicated" in line):
                b = line.split('\t')[0]
        
                if "Duplicated" in line:
                    if b not in seen:
                        seen[b] = 0
                    seen[b] += 1
                
                    name = prefix + '-' + spname + '-' + b + '-D' + str(seen[b])
                else:
                    name = prefix + '-' + spname + '-' + b + '-S1'
                dest.write(line.split('\t')[2] + '\t' + line.split('\t')[3] + '\t' + line.split('\t')[4] + '\t' + name + '\ta1000\t' + line.split('\t')[5] + '\n')
                #dest.write('\t'.join(line.split('\t')[2:6]) + '\t' + name + '\n')
    
if __name__ == "__main__":
  main(sys.argv)
