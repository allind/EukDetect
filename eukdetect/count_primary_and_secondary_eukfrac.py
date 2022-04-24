#! /usr/bin/env python
#usage: script.py
#from memory_profiler import profile
from ete3 import NCBITaxa
import argparse
import textwrap
import sys
import re

#@profile
def main(argv):

	parser = argparse.ArgumentParser(
		description=textwrap.dedent("""\
			Summarize and filter alignments by taxid.

			Required arguments are --dbfile, --inherited_markers, --taxid_link, --readcounts, --primarytab, --eukfrac, --alltab, --taxid_genelens

			"""),
		formatter_class = argparse.RawDescriptionHelpFormatter
		)

	parser.add_argument(
		"--dbfile",
		type=str,
		action="store",
		dest="dbfile",
		help= "Eukdetect database folder",
		required=True
		)


	parser.add_argument(
		"--inherited_markers",
		type=str,
		action="store",
		dest="inherited_markers",
		help= "Eukdetect database folder",
		required=True
		)

	parser.add_argument(
		"--taxid_link",
		type=str,
		action="store",
		dest="taxid_link",
		help= "Eukdetect database folder",
		required=True
		)


	parser.add_argument(
		"--readcounts",
		type=str,
		action="store",
		dest="readcounts",
		help= "Read counts and mismatches file.",
		required=True
		)

	parser.add_argument(
		"--eukfrac",
		type=str,
		action="store",
		dest="eukfrac",
		help= "Eukaryotic abundance & fraction output file",
		required=True
		)

	parser.add_argument(
		"--primarytab",
		type=str,
		action="store",
		dest="primarytab",
		help= "Table output of filtered hits.",
		required=True
		)

	parser.add_argument(
		"--alltab",
		type=str,
		action="store",
		dest="alltab",
		help= "Table output of all hits.",
		required=True
		)
	parser.add_argument(
			"--taxid_genelens",
			type = str,
			action = "store",
			dest = "taxid_genelens",
			help = "Cumulative gene length per taxid",
			required = True
			)

	files = parser.parse_args()

	#initialize NCBI taxdb
	ncbi = NCBITaxa(files.dbfile)


	#taxid genelength correspondence
	#taxid_genelen = {taxid: length}
	taxid_genelen = {line.split('\t')[0]: int(line.split('\t')[1].strip('\n')) for line in open(files.taxid_genelens)}

	#create 2 dicts for ease of lookup

	#correspondence between taxid & marker gene name
	#taxid_seqs: {taxid: [seq1, seq2]}. Save every seen taxid and which seqs
	#seq_taxids = {seq: taxid, seq:taxid} Save every seq

	taxid_seqs = {}
	seq_taxids = {}
	for line in open(files.taxid_link):
		line = line.strip('\n')
		taxid = line.split('\t')[1]
		if taxid not in taxid_seqs:
			taxid_seqs[taxid] = []
		seq = line.split('\t')[0]
		taxid_seqs[taxid].append(seq)
		seq_taxids[seq] = taxid

	#save contents of read_counts_and_mismatches file as dict per observed taxid
	#save observed genuses
	#taxid_counts: {taxid: [[marker, readcount, correct_bases, total_bases, seqlen, coverage, pid, busco]]}
	taxid_counts = {}

	counter = 0
	countfile = open(files.readcounts)
	countfile.readline()

	genuses = {}
	above_species = []
	#genuses: {genus:[taxid, taxid, taxid]}

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
		pid = float(line.split('\t')[7])
		taxid = seq_taxids[seq]

		if "Collapse" not in seq:
			busco = re.findall('-\d*at\d*-', seq)[0].strip('-')
		else:
			busco = "Collapsed"
		
		#determine genus
		lineage = ncbi.get_lineage(int(taxid))
		ranks = {value: key for (key, value) in ncbi.get_rank(lineage).items()}
		#lowest = list(ncbi.get_rank([lineage[-1]]).values())[0]
		if 'genus' in ranks and 'Collapse' not in seq and "species" in ranks: #lowest != "genus": #dont filter if it's at the genus level
			genus = ranks['genus']
			if genus not in genuses:
				genuses[genus] = []
			if taxid not in genuses[genus]:
				genuses[genus].append(taxid)
		elif "SSCollapse" not in seq: #don't add anything that's got SSCollapse in it
			above_species.append(taxid)

		#save info per sequence in seq_counts dict
		#seq_counts[seq] = [count, correct_bases, total_bases, subjlen, coverage, pid, busco]

		if taxid not in taxid_counts:
			taxid_counts[taxid] = []
			#find the genus if not a spcollapsed gene

		taxid_counts[taxid].append([seq, 
									count, 
									correct_bases, 
									total_bases, 
									subjlen, 
									coverage,
									pid,
									busco])

	if counter == 0:
		message = "Empty read count file. Likely no aligned reads in sample."
		#print(message)
		#still have to write stuff
		f = open(files.eukfrac, 'w')
		f.write(message + '\n')
		f.close()
		f = open(files.alltab, 'w')
		f.write(message + '\n')
		f = open(files.primarytab, 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()
	countfile.close()


	#done parsing read_counts_and_mismatches file

	#calculate stats for each observed taxid
	taxon_coverage = {}
	
	#taxon_coverage[taxon] = [
	#observed_markers, 
	#readcounts, 
	#total_bases, 
	#percentage_markers, 
	#marker_coverage, 
	#percent_id,
	#total_observed_marker_len,
	#buscos,
	#total_gene_length,
	#total_markers]

	seen_taxids = []
	for tax in taxid_counts:
		mc = len(taxid_counts[tax])
		counts = 0
		bases = 0
		correct = 0
		total_bases = 0
		subj_len = 0
		buscos = []
		for i in range(0, len(taxid_counts[tax])):

			busco = taxid_counts[tax][i][-1]
			if len(busco) > 1:
				buscos.append(busco)

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
		
		taxid_len = taxid_genelen[tax]
		rpkg = counts / (taxid_len/1000)

		if tax not in seen_taxids:
			seen_taxids.append(tax)

		taxon_coverage[tax] = [mc, 
								counts, 
								total_bases, 
								marker_percentage, 
								overall_coverage, 
								percent_identity, 
								subj_len,
								buscos,
								taxid_len,
								total_markers]	
	
	#create tree structure for all observed taxids

	tree = ncbi.get_topology(seen_taxids)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	tree_taxids = seen_taxids + lineage
	full_tree = ncbi.get_topology(tree_taxids, intermediate_nodes=True)
	full_taxid_lineage = [node.name for node in full_tree.traverse()]

	#full_seq_taxids: {taxid: [[specific buscos], specific count, specific + inherited count]}
	full_seq_taxids = {}
	for line in open(files.inherited_markers):
		line = line.strip('\n')
		taxid = line.split('\t')[0]
		if taxid in full_taxid_lineage:
			buscos = []
			for seq in line.split('\t')[1].split(','):
				if len(re.findall('-\d*at\d*-', seq)) > 0:
					busco = re.findall('-\d*at\d*-',seq)[0].strip('-')
					if busco not in buscos:
						buscos.append(busco)

			specific_count = len(line.split('\t')[1].split(','))
			sp_and_inherited_count = len(line.split('\t')[2].split(','))

			full_seq_taxids[taxid] = [buscos, specific_count, sp_and_inherited_count]

	#write full table
	marker_sorted = sorted(taxon_coverage.keys(), reverse = True, key = lambda x: taxon_coverage[x][3])

	dest = open(files.alltab, 'w')
	dest.write("Name\tTaxid\tRank\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\tAmount of marker length in EukDetect db\n")
	for tax in marker_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]

		if rank == "no rank":
			#parent rank
			parent = ncbi.get_lineage(tax)[-2]
			rank = [ncbi.get_rank([parent])[e] for e in ncbi.get_rank([parent])][0]			

		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		total_marker_len = taxon_coverage[tax][6]
		blen = taxid_genelen[tax]
		dest.write(name + '\t'
			+ str(tax) + '\t'
			+ rank + '\t'
			+ str(mc) + '\t' 
			+ str(counts) + '\t' 
			+ str(marker_percentage) + '%\t'
			+ str(overall_coverage) + '%\t'
			+ str(percent_identity) + '%\t'
			+ str(blen) + '\n')
	dest.close()


	#determine primary and secondary hits
	#if MRCA is at the level of genus, consider whether one should be primary or secondary by looking at buscos
	primary = {}
	secondary = {}
	genus_secondary_hits = {} #structure: {genus: [secondary_hit_taxid]}

	for g in genuses:
		if len(genuses[g]) > 1: #multiple species in same genus
			taxids = genuses[g]
			reads = [taxon_coverage[taxid][1] for taxid in taxids]
			bases = [taxon_coverage[taxid][2] for taxid in taxids]

			#if one has more reads and more bases than all others, it is primary, others are secondary
			maxreads = max(reads)
			maxbases = max(bases)
			ptaxids = []

			if (reads.count(maxreads) == 1 and bases.count(maxbases) == 1)\
			 and (reads.index(maxreads) == bases.index(maxbases)): #no ties, same ID
			 	maxtax = taxids[reads.index(maxreads)]
			 	primary[maxtax] = taxon_coverage[maxtax][0:5]
			 	ptaxids.append(maxtax)
				#ptaxids.append(taxids[reads.index(maxreads)])
				#primary[ptaxid] = taxon_coverage[ptaxid][0:5]
				#p_buscos = full_seq_taxids[ptaxid][0]
			else:
				for t in taxids: 
					if taxon_coverage[t][1] == maxreads or taxon_coverage[t][2] == maxbases:
						ptaxids.append(t)
						primary[t] = taxon_coverage[t][0:5]

			unsorted_ataxids = [t for t in taxids if t not in ptaxids]
			ataxids = sorted(unsorted_ataxids, key = lambda x: taxon_coverage[x][1], reverse = True)
			for ataxid in ataxids:

				is_secondary = False
				for ptaxid in primary:
					p_buscos = [b for b in full_seq_taxids[ptaxid][0]]
					a_buscos = taxon_coverage[ataxid][7]
					a_remain = [b for b in a_buscos if b in p_buscos]

					if len(a_remain) > 0:
						a_above = []
						for b in a_remain:

							#it may not be a hit for the other one! check first
							#check that the pid for this hit is lower
							apid = [seq[6] for seq in taxid_counts[ataxid] if seq[7] == b]
							ppid = [seq[6] for seq in taxid_counts[ptaxid] if seq[7] == b]

							if len(ppid) > 0 and apid[0] >= ppid[0]:
								a_above.append(b)
							elif len(ppid) == 0:
								a_above.append(b)
						#if a_buscos is fewer than 5, all must be correct
						#print(a_above)
						if len(a_buscos) < 5:
							if len(a_above) < len(a_buscos):
								is_secondary = True
						else:
							if len(a_above) <= len(a_buscos)/2: #change: alt hit has to be half or busco hits being above
								is_secondary = True
					else:
						is_secondary = True
				if is_secondary:
					secondary[ataxid] = taxon_coverage[ataxid][0:5] + [ptaxid]
					genus = str(g)
					if genus not in genus_secondary_hits:
						genus_secondary_hits[genus] = []
					genus_secondary_hits[genus].append(ataxid)
					#secondary_hit_reads[g].append([ataxid, taxon_coverage[ataxid][1], taxid_genelen[ataxid]])
				else:
					primary[ataxid] = taxon_coverage[ataxid][0:5]
		else: #primary
			taxid = genuses[g][0]
			primary[taxid] = taxon_coverage[taxid][0:5]

	#add anything else
	for t in above_species:
		primary[t] = taxon_coverage[t][0:5]


	primary_sorted = sorted(primary.keys(), reverse = True, key = lambda x: primary[x][3])
	#secondary_sorted = sorted(secondary.keys(), reverse=True, key=lambda x: secondary[x][3])

	filter_passing_taxids = []

	for tax in primary_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		if rank == "no rank":
			prev = ncbi.get_lineage(tax)[-1]
			prevrank = [ncbi.get_rank([prev])[e] for e in ncbi.get_rank([prev])][0]
			if prevrank == "species":
				rank = "species"
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]

		#filter
		if int(mc) >= 2 and int(counts) >= 4:
			filter_passing_taxids.append(tax)


	#close if no filter passing taxids
	if len(filter_passing_taxids) == 0:
		message = "No taxa passing filter requirements."
		#print(message)

		#still have to write stuff
		f = open(files.primarytab, 'w')
		f.write(message + '\n')
		f.close()
		f = open(files.eukfrac, 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()


	#create NCBI taxon tree of observed taxa + extend to cellular_org
	tree = ncbi.get_topology(filter_passing_taxids)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	primary_tree_taxids = [int(e) for e in filter_passing_taxids] + lineage
	primary_tree = ncbi.get_topology(primary_tree_taxids, intermediate_nodes=True)


	orphan_children = []


	#phylum class order family genus species

	taxid_lendenoms = {} #for all species, get full marker possibilities, for higher rank, get just what's specific

	#find counts of seqs for internal nodes
	relab_levels = {'species': [], 'genus': [], 'family': [], 'order': [], 'class': [], 'phylum': []}
	ordered_labels = ["phylum", "class", "order", "family", "genus", "species"]

	lineages = {}


	#pre-add secondary hits to each genus
	for g in genus_secondary_hits:
		for s in genus_secondary_hits[g]:
			for seq in taxid_counts[s]:
				#if g not in taxid_lendenoms:
				#	taxid_lendenoms[g] = 0
				#taxid_lendenoms[g] += seq[4]
				if g not in taxid_counts:
					taxid_counts[g] = []
				taxid_counts[g].append(seq)

	#calculate seqs and seqlens for each taxonomic node. this goes top-down
	for node in primary_tree.traverse():

		#get lineage name
		lin_name = ""
		currname = [ncbi.get_taxid_translator([node.name])[e] for e in ncbi.get_taxid_translator([node.name])][0]
		lineage = ncbi.get_lineage(node.name)
		names = ncbi.get_taxid_translator(ncbi.get_lineage(node.name))
		ranks = ncbi.get_rank(ncbi.get_lineage(node.name))
		ranks_rev = {ranks[e]:e for e in ranks}


		#print(lineage)
		#print(ranks)
		#print(ranks_rev)
		prev_rank = ranks[lineage[-2]]


		for i in ordered_labels:
			if i in ranks_rev:
				lin_name += i + "-" + names[ranks_rev[i]] + "|"
		lin_name = lin_name.strip('|').replace(' ', "_")
		lineages[node.name] = lin_name

		#init taxid_lendenoms and taxid_counts if does not exist
		if node.name not in taxid_lendenoms:
			taxid_lendenoms[node.name] = 0
		if node.name not in taxid_counts:
			taxid_counts[node.name] = []



		rank = [ncbi.get_rank([node.name])[e] for e in ncbi.get_rank([node.name])][0]



		#if node.is_leaf() == False and rank != "species":
		#if rank != "species":	
		if rank in relab_levels:
			relab_levels[rank].append(node.name)

		if rank != "species" and prev_rank != "species": #if not a species or a strain, add individuals
			#add indiv seqs
			for seq in taxid_counts[node.name]:
				taxid_lendenoms[node.name] += seq[4]

		if (rank == "species" or prev_rank == "species") and node.name in taxid_genelen:
			taxid_lendenoms[node.name] += taxid_genelen[node.name]	

		for desc in node.iter_descendants():
			if desc.name in taxid_counts:
				descrank = [ncbi.get_rank([desc.name])[e] for e in ncbi.get_rank([desc.name])][0]
				
				dlineage = ncbi.get_lineage(node.name)
				dnames = ncbi.get_taxid_translator(ncbi.get_lineage(node.name))
				dranks = ncbi.get_rank(ncbi.get_lineage(node.name))
				d_prev_rank = dranks[dlineage[-2]]
			

				if descrank == "species" or d_prev_rank == "species":
					taxid_lendenoms[node.name] += taxid_genelen[desc.name] #if sp add full markers
				else:
					for seq in taxid_counts[desc.name]:
						taxid_lendenoms[node.name] += seq[4] #if not sp add submarkers

				for seq in taxid_counts[desc.name]:
					if seq not in taxid_counts[node.name]:
						taxid_counts[node.name].append(seq)
		#elif node.is_leaf() == False and rank == "species":
			#stuff
	#		x = 0
	#	else: #has to be a strain?

			#add full seq since is species
			#if node.name in taxid_genelen: #avoids case where there is a strain without dedicated taxid

	#		taxid_lendenoms[node.name] += taxid_genelen[node.name]
	#		relab_levels['species'].append(node.name)
	#		if node.name not in taxid_counts:
	#			orphan_children.append(node.name)


	#determine if all hits have all levels


	levels_to_remove = []
	for tax in filter_passing_taxids:
		lin = lineages[tax]
		groups = [l.split('-')[0] for l in lin.split('|')]
		levels = ordered_labels[0:len(groups)]

		if levels != groups:
			#for g in groups:
			for l in levels:
				if l not in groups:
					if l not in levels_to_remove:
						levels_to_remove.append(l)


	for l in levels_to_remove:
		relab_levels.pop(l)



	#calculate relabs for each level
	relabs = {} #relabs[taxid] = [reads, amt_marker_sequence, rpks, eukfrac]
	for group in relab_levels:
		sum_rpks = 0
		for tax in relab_levels[group]:
			reads = 0
			for seq in taxid_counts[tax]:
				reads += seq[1]
			amt_marker_sequence = taxid_lendenoms[tax]
			rpks = reads / (amt_marker_sequence / 1000)
			sum_rpks += rpks
			relabs[tax] = [reads, amt_marker_sequence, rpks]
		for tax in relab_levels[group]:
			eukfrac = (relabs[tax][2] / sum_rpks) * 100
			relabs[tax].append(eukfrac)



	dest = open(files.primarytab, 'w')
	dest.write("Name\tRank\tLineage\tTaxid\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	for tax in filter_passing_taxids:
		if tax in lineages:
			lin = lineages[tax]
		else:
			lin = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]

		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]


		if rank == "no rank":
			#parent rank
			parent = ncbi.get_lineage(tax)[-2]
			prevrank = [ncbi.get_rank([parent])[e] for e in ncbi.get_rank([parent])][0]			
			if prevrank == "species":
				rank = "species"

		#lin = lineages[tax]
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]

		#filter

		dest.write(name + '\t'
			+ rank + '\t'
			+ lin + '\t'
			+ str(tax) + '\t'
			+ str(mc) + '\t' 
			+ str(counts) + '\t' 
			+ str(marker_percentage) + '%\t'
			+ str(overall_coverage) + '%\t'
			+ str(percent_identity) + '%\n')
	dest.close()



	#table with relative abundance for all levels
	dest = open(files.eukfrac, 'w')
	dest.write("Lineage\tRank\tName\tTaxID\tRPKS\tEuk_fraction\tReads\tAmt_marker_sequence\n")

	for node in primary_tree.traverse("preorder"):
		rank = [ncbi.get_rank([node.name])[e] for e in ncbi.get_rank([node.name])][0]
		name = [ncbi.get_taxid_translator([node.name])[e] for e in ncbi.get_taxid_translator([node.name])][0]
		if node.name in lineages:
			lin = lineages[node.name]
		else:
			lin = [ncbi.get_rank([node.name])[e] for e in ncbi.get_rank([node.name])][0]

		if rank == "no rank" and node.is_leaf():
			continue #is strain, have already printed species at this point
			#rank = "species"
		if node.name in relabs:
			rpks = round(relabs[node.name][2], 4)
			eukfrac = round(relabs[node.name][3], 4)
			reads = relabs[node.name][0]
			markerseq = relabs[node.name][1]
			dest.write(lin + '\t'
				+ rank + '\t'
				+ name + '\t' 
				+ node.name + '\t'
				+ str(rpks) + '\t'
				+ str(eukfrac) + '\t'
				+ str(reads) + '\t' 
				+ str(markerseq) + '\n'
				)
	dest.close()



if __name__ == "__main__":
  main(sys.argv)
