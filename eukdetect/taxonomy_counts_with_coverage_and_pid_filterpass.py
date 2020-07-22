#! /usr/bin/env python
#usage: script.py [taxid_link] [taxa.sqlite] [inherited_markers] [idxstats] [taxonomy output] [table output]
from ete3 import NCBITaxa
import sys

def main(argv):
		
	#read in taxonomy info for each BUSCO
	species_taxids = [] #species_taxids[marker_id] = taxid

	for line in open(sys.argv[1]):
		tax = line.split('\t')[1].strip('\n')
		if tax not in species_taxids:
			species_taxids.append(tax)
			
	#initialize NCBI taxdb
	ncbi = NCBITaxa(sys.argv[2])

	#create 2 dicts for ease of lookup
	#taxid_seqs: {taxid: [seq1, seq2]}. Save every seen taxid and which seqs
	#seq_taxids = {seq: taxid, seq:taxid} Save every seq
	taxid_seqs = {}
	seq_taxids = {}
	for line in open(sys.argv[1]):
		line = line.strip('\n')
		taxid = line.split('\t')[1]
		if taxid not in taxid_seqs:
			taxid_seqs[taxid] = []
		seq = line.split('\t')[0]
		taxid_seqs[taxid].append(seq)
		seq_taxids[seq] = taxid

	#iterate over idxstats file and save counts
	#seq_counts[seq] = [readcount, correct_bases, total_bases, seqlen, coverage]
	seq_counts = {}
	
	counter = 0
	countfile = open(sys.argv[4])
	countfile.readline()
	for line in countfile:
		counter += 1
		line = line.strip('\n')
		seq = line.split('\t')[0]
		count = int(line.split('\t')[1])
		correct_bases = int(line.split('\t')[2])
		incorrect_bases = int(line.split('\t')[3])
		total_bases = int(line.split('\t')[4])
		subjlen = int(line.split('\t')[5])
		coverage = float(line.split('\t')[6])
		seq_counts[seq] = [count, correct_bases, total_bases, subjlen, coverage]
		taxid = seq_taxids[seq]
		#if taxid not in seen_taxids:
		#	seen_taxids.append(int(taxid))

	if counter == 0:
		message = "Empty read count file. Likely no aligned reads in sample."
		print(message)
		#still have to write stuff
		f = open(sys.argv[5], 'w')
		f.write(message + '\n')
		f.close()
		f = open(sys.argv[6], 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()
	#done parsing idxstats file
	


	full_seq_taxids = {line.split('\t')[0]: [line.split('\t')[1].split(','), line.split('\t')[-1].strip('\n').split(',')] for line in open(sys.argv[3])}
	#full_seq_taxids: {taxid: [[specific buscos], [specific + inherited buscos]]}
	#determine seq counts 
	
	#taxid_counts: {taxid: [[marker, readcount, correct_bases, total_bases, seqlen, coverage]]}
	taxid_counts = {}
	for seq in seq_counts:
		taxid = seq_taxids[seq]
		if taxid not in taxid_counts:
			taxid_counts[taxid] = []
		taxid_counts[taxid].append([seq, int(seq_counts[seq][0]), int(seq_counts[seq][1]), seq_counts[seq][2], seq_counts[seq][3], seq_counts[seq][4]])

	#write just observed taxid seqs
	taxon_coverage = {}
	#taxon_coverage[taxon] = [observed_markers, readcounts, total_bases, percentage_markers, marker_coverage, percent_id ]
	#dest = open(sys.argv[6], 'w')
	#dest.write("Name\tNCBI_Rank\tTaxID\tObserved_markers\tRead_counts\tPercent_observed_markers\tMarker_coverage\tPercent_identity\n")
	seen_taxids = []
	for tax in taxid_counts:
		mc = len(taxid_counts[tax])
		counts = 0
		bases = 0
		correct = 0
		total_bases = 0
		subj_len = 0
		for i in range(0, len(taxid_counts[tax])):
			counts += taxid_counts[tax][i][1]
			bases += taxid_counts[tax][i][3]
			correct += taxid_counts[tax][i][2]
			total_bases += taxid_counts[tax][i][3]
			subj_len += taxid_counts[tax][i][4]
		percent_identity = round((correct / total_bases) * 100, 2)
		overall_coverage = round((total_bases / subj_len ) * 100, 2)
		total_markers = len(taxid_seqs[tax])
		marker_percentage = round( mc / total_markers * 100, 2)
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		#rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		#here is where the filter happens:
		if mc >= 2 and counts >=4:
			taxon_coverage[tax] = [mc, counts, total_bases, marker_percentage, overall_coverage, percent_identity]
			seen_taxids.append(tax)


	if len(seen_taxids) == 0:
		message = "No taxa passing filter requirements."
		print(message)
		#still have to write stuff
		f = open(sys.argv[5], 'w')
		f.write(message + '\n')
		f.close()
		f = open(sys.argv[6], 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()

	#dest.close()
	dest = open(sys.argv[6], 'w')
	dest.write("Name\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	marker_sorted = sorted(taxon_coverage.keys(), reverse = True, key = lambda x: taxon_coverage[x][3])

	for tax in marker_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		dest.write(name + '\t'
			+ str(mc) + '\t' 
			+ str(counts) + '\t' 
			+ str(marker_percentage) + '%\t'
			+ str(overall_coverage) + '%\t'
			+ str(percent_identity) + '%\n')

	#create NCBI taxon tree of observed taxa + extend to cellular_org
	tree = ncbi.get_topology(seen_taxids)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	full_taxids = seen_taxids + lineage
	full_tree = ncbi.get_topology(full_taxids, intermediate_nodes=True)
	orphan_children = []

	#find counts of seqs for internal nodes
	for node in full_tree.traverse():
		if node.is_leaf() == False:
			if node.name not in taxid_counts:
				taxid_counts[node.name] = []
			for desc in node.iter_descendants():
				if desc.name in taxid_counts:
					for seq in taxid_counts[desc.name]:
						if seq not in taxid_counts[node.name]:
							taxid_counts[node.name].append(seq)
		else:
			if node.name not in taxid_counts:
				orphan_children.append(node.name)
	


	#print the tree	
	level_counts = []
	currspaces = 0
	currparent = ''
	seen_parents = {}
	dest = open(sys.argv[5], 'w')
	dest.write("Markers_Obs\tTotal_Markers\tPercent_Makers_Obs\tPercent_ID\tMarker_read_count\tRank\tName\n")
	for node in full_tree.traverse("preorder"):
		if node.name not in orphan_children:
			rank = [ncbi.get_rank([node.name])[e] for e in ncbi.get_rank([node.name])][0]
			name = [ncbi.get_taxid_translator([node.name])[e] for e in ncbi.get_taxid_translator([node.name])][0]
			if node.is_root():
				currspaces = 0
			else:
				if currparent == '':
					currparent = node.up.name
					currspaces += 4
				else:
					if currparent != node.up.name:
						currparent = node.up.name
						if currparent in seen_parents:
							currspaces = seen_parents[currparent]
						else:
							currspaces += 4
							seen_parents[currparent] = currspaces
			if node.name in taxon_coverage:
				pid = str(taxon_coverage[node.name][5]) + '%'
			else:
				pid = "NA"
			#total_buscos
			buscos = len(taxid_counts[node.name])
			seqs = sum([b[1] for b in taxid_counts[node.name]])
			total_buscos = len(full_seq_taxids[node.name][1])
			percent = round((buscos/total_buscos)*100,2)
			dest.write(str(buscos) + '\t' 
				+ str(total_buscos) + "\t" 
				+ str(percent)  + '%\t' 
				+ str(pid) + '\t'
				+ str(seqs) + '\t' 
				+ rank + '\t' 
				+ ' ' * currspaces + name + '\n')
	dest.close()

if __name__ == "__main__":
  main(sys.argv)
