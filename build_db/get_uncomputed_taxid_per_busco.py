#! /usr/bin/env python
#usage: script.py species_taxids fasta collapsed_ids
from Bio import SeqIO
from ete3 import NCBITaxa
import sys, re, argparse
def main(argv):
	parser = argparse.ArgumentParser(description="Arguments for computing taxids")
	parser.add_argument("--taxdb",required=True,  type=str, help="Path to taxa.sqlite (taxa.sqlite.traverse.pkl assumed to be in same dir)")
	parser.add_argument("--speciestax", required=True, type=str, help="Tab delimited file of species name (as encoded in busco header) and taxonomy ID")
	parser.add_argument("--fasta", required=True, type=str, help="Final busco fasta")
	parser.add_argument("--collapsed_ids", required=True, type=str, help="Output from rename_cdhit_collapsed.py script - named something like *_collapsed_seqnames.txt")
	
	args = parser.parse_args()
	taxdbpath = args.taxdb
	speciestax = args.speciestax
	fasta = args.fasta
	collapsed = args.collapsed_ids

	all_taxids = {}
	notfound = []
	ncbi = NCBITaxa(dbfile = taxdbpath)
	sp_taxids = {line.split('\t')[0]: line.split('\t')[1].strip('\n') for line in open(speciestax)}
	#known_taxids = {line.split('\t')[0]: line.split('\t')[1].strip('\n') for line in open("previous_busco_taxid_link.txt")}

	#all_taxids: {seq: taxid, seq:taxid}
	for seq in SeqIO.parse(fasta, 'fasta'):
		sp = re.split('-\d*at\d*-', '-'.join(seq.id.split('-')[1:]))[0]
		#if seq.id not in known_taxids:
		if sp not in sp_taxids:
			notfound.append(seq.id)
		else:
			all_taxids[seq.id] = sp_taxids[sp]
		#else:
		#	all_taxids[seq.id] = known_taxids[seq.id]
	counter = 0

	for line in open(collapsed):
		counter += 1
		print(counter)
		line = line.strip('\n')
		seq = line.split('\t')[0]
		other_buscos = line.split('\t')[3:-1]
		other_species = []
		for sp in other_buscos:
			#new = sp.split("_EO")[0].replace("_", " ")
			new = re.split('-\d*at\d*-', '-'.join(sp.split('-')[1:]))[0]
			
			if " sp " in new:
				new = new.replace(" sp ", " sp. ")
			other_species.append(new)
		#taxids = ncbi.get_name_translator(other_species)
		taxids = [sp_taxids[sp] for sp in other_species]
		#print(other_species)
		#print(taxids)
		#for sp in other_species:
		#	if sp not in taxids.keys():
		#		if len(sp.split(' ')) > 2:
		#			new = ' '.join(sp.split(' ')[0:2])
		#			newtaxid = ncbi.get_name_translator([new])
		#			if len(newtaxid.keys()) == 0:
		#				newtaxid = ncbi.get_name_translator([new.split(' ')[0]])
		#			taxids.update(newtaxid)
		#	#find mrca
		#print(taxids)
		#taxid_nums = [taxids[taxid][0] for taxid in taxids]
		#print(taxid_nums)
		tree = ncbi.get_topology(taxids)
		#print(other_species)
		#print(tree.get_ascii(attributes=["sci_name", "rank"]))
		#print(taxid_nums)
		#print("Root", tree.get_tree_root().name)
		root = tree.get_tree_root().name
		all_taxids[seq] = root
		#print(root)
		#print(seq + '\t' + root)
	for t in all_taxids:
		print(t + '\t' + all_taxids[t])
if __name__ == "__main__":
  main(sys.argv)
