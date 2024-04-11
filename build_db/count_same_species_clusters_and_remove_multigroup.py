#! /usr/bin/env python
import sys, re
def main(argv):
    clusters = {}
    cluster = ""
    cluster_mains = {}
    for line in open(sys.argv[1]):
        line = line.strip('\n')
        if ">" in line[0]:
            cluster = line.strip('>')
            clusters[cluster] = []
        else:
            clusters[cluster].append(line.split('>')[1].split('..')[0])

        if "... *" in line:
            cluster_mains[cluster] = line.split('>')[1].split('..')[0]
    same_species = 0
    same_genus = 0
    other = 0
    remove = []
    for cluster in clusters:
        species = []
        genus = []
        groups = []
        toprint = True
        collapsed = False
        symbiont = False
        for seq in clusters[cluster]:
            #sp = seq.split('_EO')[0]
            sp = re.split('-\d*at\d*-', '-'.join(seq.split('-')[1:]))[0]
            species.append(sp)
            g = sp.split("_")[0]
            genus.append(g)
            group = seq.split('-')[0]
            groups.append(group)
            if "Collapse" in seq:
                collapsed = True
                toprint = False
            if "symbiont" in sp:
                symbiont = True
        unique_groups = list(set(groups))
        unique_species = (list(set(species)))
        unique_genus = list(set(genus))
        name = ""
        if len(unique_species) == 1 and len(unique_groups) == 1:
            same_species += 1
            name = "same_species"
        elif len(unique_genus) == 1 and len(unique_groups) == 1:
            same_genus += 1
            name = "same_genus"
        elif len(unique_groups) > 1 and collapsed == False:
            toprint = False
            remove.append(cluster)
        else:
            other += 1
            name = "other"
        if toprint == True:
            print(cluster_mains[cluster] + '\t' + cluster + '\t' + "\t".join(clusters[cluster]) + '\t' + name)
    #print("Same species:", same_species)
    #print("Same genus:", same_genus)
    #print("Other:", other)
    prefix=sys.argv[1].split('.')[0]
    dest = open(prefix + "_removed_for_multigroup.txt", 'w')
    for r in remove:
        #dest.write(cluster_mains[r] + '\t' + r + '\t' + "\t".join(clusters[r]) + '\n')
        dest.write('\n'.join(clusters[r]) + '\n')
    dest.close()

if __name__ == "__main__":
  main(sys.argv)
