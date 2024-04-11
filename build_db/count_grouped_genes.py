#! /usr/bin/env python
from Bio import SeqIO
import sys,re
def main(argv):


    gene_ids = {} #gene_ids[species] = [num_single_buscos,num_duplicated_buscos,total_num_genes,num_grouped_buscos,{species_grouped: count of genes, species_grouped: count of genes}]
    #0 is num_single, 1 is num_dup, 2 is total_num_gnees, 3 is collapsed_buscos, 4 is dict with species identifier of grouped genes and their counts
    for seq in SeqIO.parse(sys.argv[1], 'fasta'):
        sp = re.split('-\d*at\d*-', '-'.join(seq.id.split('-')[1:]))[0]
        status = seq.id.split('-')[-1]
        if sp not in gene_ids:
            gene_ids[sp] = [0,0,0,0, {}]
        gene_ids[sp][2] += 1
        if "D" in status:
            gene_ids[sp][1] += 1
        else:
            gene_ids[sp][0] += 1

    for line in open(sys.argv[2]):
        line = line.strip('\n')
        buscos = line.split('\t')[2:-1]
        species = [re.split('-\d*at\d*-', '-'.join(e.split('-')[1:]))[0] for e in buscos]
        uniq_species = list(set(species))
        spcounts = {s: species.count(s) for s in uniq_species}
        for s in species:
            gene_ids[s][3] += 1
            for os in spcounts:
                if os != s:
                    if os not in gene_ids[s][4]:
                        gene_ids[s][4][os] = 0
                    gene_ids[s][4][os] += spcounts[os]
                else:
                    if os not in gene_ids[s][4]:
                        gene_ids[s][4][os] = 0
                    gene_ids[s][4][os] += spcounts[os] - 1

        



    print("Species\tSingle_copy_buscos\tDuplicated_buscos\tTotal_buscos\tCollapsed_buscos\tSpecies_collapsed")
    for g in gene_ids:
        print(g + '\t' + "\t".join([str(e) for e in gene_ids[g][0:4]]) + '\t' + ",".join([e + ":" + str(gene_ids[g][4][e]) for e in gene_ids[g][4]]))

            
    
if __name__ == "__main__":
  main(sys.argv)
