#! /usr/bin/env python
#usage: script.py [busco_taxid_link] [taxa.sqlite]
from ete3 import NCBITaxa
#from memory_profiler import profile
import sys

#@profile
def main(argv):
		
#read in taxonomy info for each BUSCO
	species_taxids = [] #species_taxids[marker_id] = taxid

	for line in open(sys.argv[1]): #parse busco_taxid_link
		tax = line.split('\t')[1].strip('\n')
		if tax not in species_taxids:
			species_taxids.append(tax)


	#initialize NCBI taxdb
	ncbi = NCBITaxa(dbfile = sys.argv[2])
	
	#find the lineage to treet root
	tree = ncbi.get_topology(species_taxids)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)

	#propagate marker tree back to cellular_organisms
	full_taxids = species_taxids + lineage
	euk_tree = ncbi.get_topology(full_taxids, intermediate_nodes=True)

	#create 2 dicts for ease of lookup
	#taxid_seqs: {taxid: [seq1, seq2]}. Save every seen taxid and which seqs
	
	###seq_taxids = {seq: taxid, seq:taxid} Save every seq

	taxid_seqs = {}
	#seq_taxids = {}
	for line in open(sys.argv[1]):#parse busco_taxid_link again
		line = line.strip('\n')
		taxid = line.split('\t')[1]
		if taxid not in taxid_seqs:
			taxid_seqs[taxid] = []
		seq = line.split('\t')[0]
		taxid_seqs[taxid].append(seq)
		#seq_taxids[seq] = taxid

	#traverse tree to figure out full possible BUSCOs
	full_seq_taxids = {}
	#full_seq_taxids: {taxid: [[specific buscos], [specific + inherited buscos]]}

	dest = open("specific_and_inherited_markers_per_taxid.txt", 'w')
	dest.write("TaxID\tSpecific_markers\tDescendant_markers\n")
	#traverse full euk tree bottom up to get possible seqs at each node

	#descendants: {taxid: [list of descendants], taxid:[list of descendants]}
	descendants = {}

	for node in euk_tree.traverse():
		#realname= [ncbi.get_taxid_translator([node.name])[e] for e in ncbi.get_taxid_translator([node.name])][0]
		# if node.is_leaf():

		# 	if node.name in taxid_seqs:
		# 		#specific = inherited
		# 		full_seq_taxids[node.name] = [taxid_seqs[node.name], taxid_seqs[node.name]]
		# 	else:
		# 		#there is nothing specific or inherited
		# 		full_seq_taxids[node.name] = [[],[]] 
		# else:
		# 	#internal node
		# 	sp_seq = []
		# 	#check if node has assigned seqs
		# 	if node.name in taxid_seqs:
		# 		sp_seq = taxid_seqs[node.name]
			
		# 	full_seq_taxids[node.name] = [sp_seq]
			
			
		# dest.write(node.name + '\t' + ",".join(full_seq_taxids[node.name][0]) + '\n')
		if node.is_leaf():
			if node.name in taxid_seqs:
				#specific = inherited
				full_seq_taxids[node.name] = [taxid_seqs[node.name], taxid_seqs[node.name]]
			else:
				#there is nothing specific or inherited
				full_seq_taxids[node.name] = [[],[]] 
		else:
			#internal node
			sp_seq = []
			#check if node has assigned seqs
			if node.name in taxid_seqs:
				sp_seq = taxid_seqs[node.name]
			
			#find inherited seqs
			descent_seqs = []
			#iterate over descendants and save specific taxids
			for desc in node.iter_descendants():
				#print(descent_seqs)
				if desc.name in taxid_seqs:
					sq = taxid_seqs[desc.name]
					descent_seqs = descent_seqs + sq
					#for s in sq:
					#	if s not in descent_seqs:
					#		descent_seqs.append(s)
				else:#this means the descendant seq is not in taxid_seqs
					if desc.name not in full_seq_taxids:
						full_seq_taxids[node.name] = [[],[]]
			#update dict with internal node seqs
			full_seq_taxids[node.name] = [sp_seq, descent_seqs]
		

		dest.write(node.name + '\t' + ",".join(full_seq_taxids[node.name][0]) + '\t' + ",".join(full_seq_taxids[node.name][1]) + '\n')
		#reset inherited taxids for memory reduction
		full_seq_taxids[node.name][1] = []
		
	dest.close() 
		

if __name__ == "__main__":
  main(sys.argv)
