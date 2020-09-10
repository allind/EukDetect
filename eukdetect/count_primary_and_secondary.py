#! /usr/bin/env python
#usage: script.py [4read_counts_and_mismatches] [5taxonomy output] [6table output] [7all hits]
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

			Required arguments are --dbfile, --readcounts, --primarytax, --primarytab, --alltab

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
		"--primarytax",
		type=str,
		action="store",
		dest="primarytax",
		help= "Taxonomy output of filtered hits",
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

	files = parser.parse_args()

	#initialize NCBI taxdb
	ncbi = NCBITaxa(files.dbfile)

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
			busco = ""
		
		#determine genus
		lineage = ncbi.get_lineage(int(taxid))
		ranks = {value: key for (key, value) in ncbi.get_rank(lineage).items()}
		lowest = list(ncbi.get_rank([lineage[-1]]).values())[0]
		if 'genus' in ranks and 'Collapse' not in seq and lowest != "genus": #dont filter if it's at the genus level
			genus = ranks['genus']
			if genus not in genuses:
				genuses[genus] = []
			if taxid not in genuses[genus]:
				genuses[genus].append(taxid)

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
		f = open(files.primarytax, 'w')
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
	
	#taxon_coverage[taxon] = [observed_markers, 
	#readcounts, 
	#total_bases, 
	#percentage_markers, 
	#marker_coverage, 
	#percent_id, 
	#buscos]

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

		if tax not in seen_taxids:
			seen_taxids.append(tax)

		taxon_coverage[tax] = [mc, 
								counts, 
								total_bases, 
								marker_percentage, 
								overall_coverage, 
								percent_identity, 
								buscos]

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

	#determine primary and secondary hits
	#if MRCA is at the level of genus, consider whether one should be primary or secondary by looking at buscos
	primary = {}
	secondary = {}
	for g in genuses:
		if len(genuses[g]) > 1: #multiple species in same genus
			taxids = genuses[g]
			reads = [taxon_coverage[taxid][1] for taxid in taxids]
			bases = [taxon_coverage[taxid][2] for taxid in taxids]

			#if one has more reads and more bases than all others, it is primary, others are secondary
			maxreads = max(reads)
			maxbases = max(bases)

			if (reads.count(maxreads) == 1 and bases.count(maxbases) == 1) and (reads.index(maxreads) == bases.index(maxbases)): #no ties, same ID
				ptaxid = taxids[reads.index(maxreads)]
				primary[ptaxid] = taxon_coverage[ptaxid][0:5]
				p_buscos = full_seq_taxids[ptaxid][0]
				ataxids = [t for t in taxids if t != ptaxid]

				for ataxid in ataxids:

					a_buscos = taxon_coverage[ataxid][-1]
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

						if len(a_above) >= len(a_remain)/2: #if half or more of the hits were above the main hit in PID
							primary[ataxid] = taxon_coverage[ataxid][0:5] #cant distinguish

						else:
							secondary[ataxid] = taxon_coverage[ataxid][0:5] + [ptaxid]
					else:
						secondary[ataxid] = taxon_coverage[ataxid][0:5]

#TODO: figure out exceptions!
			else: #there are ties or max reads and max bases are different
				#if reads.index(maxreads) != bases.index(maxbases): #max reads and max bases differ
				#	ptaxids = [t for t in taxids if taxon_coverage[t][1] == maxreads or taxon_coverage[t][2] == maxbases]
				#elif reads.count(maxreads) != 1 or bases.count(maxbases) != 1: #there are ties
				ptaxids = [t for t in taxids if taxon_coverage[t][1] == maxreads or taxon_coverage[t][2] == maxbases]

				for p in ptaxids:
					primary[p] = taxon_coverage[p][0:5]


				ataxids = [t for t in taxids if t not in ptaxids]
				a_nested = []

				if len(ataxids) > 0:

					for p in ptaxids:
						p_buscos = full_seq_taxids[p][0]

						for a in ataxids:
							a_buscos = taxon_coverage[a][-1]
							a_remain = [b for b in a_buscos if b in p_buscos]

							if len(a_remain) > 0:
								a_above = []
								for b in a_remain:
									apid = [seq[6] for seq in taxid_counts[a] if seq[7] == b]
									ppid = [seq[6] for seq in taxid_counts[p] if seq[7] == b]

									if len(ppid) > 0 and apid[0] >= ppid[0]:
										a_above.append(b)
									elif len(ppid) == 0:
										a_above.append(b)

								if len(a_above) < len(a_remain)/2: #if not all all of the hits were above the main hit in PID
									#is nested
									a_nested.append(a)
									break
							else:
								a_nested.append(a)
								break

					a_passed = [a for a in ataxids if a not in a_nested]
					for a in a_passed:
						primary[a] = taxon_coverage[a][0:5]
					for a in a_nested:
						secondary[a] = taxon_coverage[a][0:5]



		else: #primary
			taxid = genuses[g][0]
			primary[taxid] = taxon_coverage[taxid][0:5]


	dest = open(files.alltab, 'w')
	dest.write("Name\tTaxid\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	marker_sorted = sorted(taxon_coverage.keys(), reverse = True, key = lambda x: taxon_coverage[x][3])
	#TODO: implement filters

	primary_sorted = sorted(primary.keys(), reverse = True, key = lambda x: primary[x][3])
	#secondary_sorted = sorted(secondary.keys(), reverse=True, key=lambda x: secondary[x][3])
	filter_passing_taxids = []

	for tax in primary_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		#filter
		if int(mc) >= 2 and int(counts) >= 4:
			filter_passing_taxids.append(tax)
			dest.write(name + '\t'
				+ str(tax) + '\t'
				+ str(mc) + '\t' 
				+ str(counts) + '\t' 
				+ str(marker_percentage) + '%\t'
				+ str(overall_coverage) + '%\t'
				+ str(percent_identity) + '%\n')
	dest.close()

	#close if no filter passing taxids
	if len(filter_passing_taxids) == 0:
		message = "No taxa passing filter requirements."
		#print(message)

		#still have to write stuff
		f = open(files.primarytab, 'w')
		f.write(message + '\n')
		f.close()
		f = open(files.primarytax, 'w')
		f.write(message + '\n')
		f.close()
		sys.exit()


	dest = open(files.primarytab, 'w')
	dest.write("Name\tTaxid\tObserved_markers\tRead_counts\tPercent_observed_markers\tTotal_marker_coverage\tPercent_identity\n")
	for tax in marker_sorted:
		rank = [ncbi.get_rank([tax])[e] for e in ncbi.get_rank([tax])][0]
		name = [ncbi.get_taxid_translator([tax])[e] for e in ncbi.get_taxid_translator([tax])][0]
		mc = taxon_coverage[tax][0]
		counts = taxon_coverage[tax][1]
		marker_percentage = taxon_coverage[tax][3]
		overall_coverage = taxon_coverage[tax][4]
		percent_identity = taxon_coverage[tax][5]
		dest.write(name + '\t'
			+ str(tax) + '\t'
			+ str(mc) + '\t' 
			+ str(counts) + '\t' 
			+ str(marker_percentage) + '%\t'
			+ str(overall_coverage) + '%\t'
			+ str(percent_identity) + '%\n')
	dest.close()


	#create NCBI taxon tree of observed taxa + extend to cellular_org
	tree = ncbi.get_topology(filter_passing_taxids)
	tree_root = tree.get_tree_root().name
	lineage = ncbi.get_lineage(tree_root)
	primary_tree_taxids = [int(e) for e in filter_passing_taxids] + lineage
	primary_tree = ncbi.get_topology(primary_tree_taxids, intermediate_nodes=True)
	#write the tree	structure to file

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


	#create new tree of filter passing hits


	level_counts = []
	currspaces = 0
	currparent = ''
	seen_parents = {}

	dest = open(files.primarytax, 'w')
	dest.write("Markers_Obs\tTotal_Markers\tPercent_Makers_Obs\tPercent_ID\tMarker_read_count\tRank\tName\n")
	for node in primary_tree.traverse("preorder"):

		if node.name not in orphan_children and node.name in full_seq_taxids:
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
			buscos = len(taxid_counts[str(node.name)])
			seqs = sum([b[1] for b in taxid_counts[node.name]])
			total_buscos = full_seq_taxids[node.name][2]
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
