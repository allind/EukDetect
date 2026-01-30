#! /usr/bin/env python
from ete3 import NCBITaxa
import sys
import argparse
import textwrap
import sys
import logging
import os


def main(argv):

	parser = argparse.ArgumentParser(
		description=textwrap.dedent("""\
			Combine eukfracs across multiple samples.

			Input is a comma separated list of sample names, the eukdetect output directory, and the eukdetect database location.

			Output is a concatenated file with RPKS and eukfrac.

			IMPORTANT NOTE: RPKS is not comparable between samples as it is not scaled by library size. It is included here to provide context about the absolute amount of eukaryotic sequence present in a library.

			"""),
		formatter_class = argparse.RawDescriptionHelpFormatter
		)

	parser.add_argument(
		"--samples",
		type=str,
		action="store",
		dest="samples",
		help= "Comma-separated list of samples to combine",
		required=True
		)


	parser.add_argument(
		"--sampledir",
		type=str,
		action="store",
		dest="sampledir",
		help="Path to eukdetect output directory",
		required=True
		)

	parser.add_argument(
		"--output",
		type=str,
		action="store",
		dest="output",
		help= "Output file",
		required=True
		)

	parser.add_argument(
		"--eukdb",
		type=str,
		action="store",
		dest="eukdb",
		help="Eukdetect database path",
		required=True
		)

	files = parser.parse_args()

	samples = files.samples.split(',')
	samples_with_hits = []

	ncbi = NCBITaxa(files.eukdb + '/taxa.sqlite')

	sample_taxids = {}



	species = {}
	lineage_taxid_link = {}
	all_taxids = []
	sample_taxids = {}

	for sample in samples:
		e = files.sampledir + '/' + sample + "_filtered_hits_eukfrac.txt"


		if not os.path.isfile(e):
			logging.error("File does not exist: ", e)
			exit(1)
		else:
			efile = open(e)
			efile.readline()
			keep = False
			if sample not in sample_taxids:
				sample_taxids[sample] = []
			for line in efile:
				keep = True

				line = line.strip('\n')
				sp = "\t".join(line.split('\t')[0:4])
				taxid = line.split('\t')[3]
				rpks = line.split('\t')[4]
				eukfrac = line.split('\t')[5]
				reads = line.split('\t')[6]
				markers = line.split('\t')[7]
				sample_taxids[sample].append(taxid)
				if taxid not in all_taxids:
					all_taxids.append(taxid)

				if taxid not in lineage_taxid_link:
					lineage_taxid_link[taxid] = sp

				if sp not in species:
					species[sp] = {}
				species[sp][sample] = [rpks, eukfrac, reads, markers]

			if keep == False:
				sample_taxids.pop(sample)
			else:
				samples_with_hits.append(sample)


	dest = open(files.output, 'w')
	
	towrite = "Lineage\tRank\tName\tTaxID\t"
	for sample in samples_with_hits:
		towrite += sample + "_RPKS\t" + sample + "_eukfrac\t" + sample + "_reads\t" + sample + "_amt_marker_sequence\t"
	towrite = towrite.strip('\t')
	towrite += "\n"
	dest.write(towrite)


	all_tree = ncbi.get_topology(all_taxids)

	#sample_trees = {sample: ncbi.get_topology(sample_taxids[sample]) for sample in sample_taxids}

	for node in all_tree.traverse("preorder"):

		if node.name in lineage_taxid_link: #was actually calculated somewhere

			sp = lineage_taxid_link[node.name]
			seen = list(species[sp].keys())
			not_seen = [s for s in samples_with_hits if s not in seen]

			for s in not_seen:
				has_descendants = False
				#figure out if its descendants are seen
				s_taxids = sample_taxids[s]
				for desc in node.iter_descendants():
					if desc.name in s_taxids:
						#print(desc.name)
						has_descendants = True


				if has_descendants:
					species[sp][s] = ["Not calculated","Not calculated", "Not calculated", "Not calculated"]
				else:
					species[sp][s] = ["0","0", "0", "NA"]

			towrite = sp + '\t'
			for sample in samples_with_hits:
				towrite += "\t".join(species[sp][sample]) + '\t'
			towrite = towrite.strip('\t')
			towrite += "\n"
			dest.write(towrite)
	dest.close()

	#for sp in species:

	#	seen = list(species[sp].keys())

	#	not_seen = [s for s in samples if s not in seen]

	#	for s in not_seen:
	#		species[sp][s] = ["0","0", "0", "NA"]

	#	towrite = sp + '\t'
	#	for sample in samples:
	#		towrite += "\t".join(species[sp][sample])

	#	towrite += "\n"

	#	dest.write(towrite)

	#dest.close()





if __name__ == "__main__":
  main(sys.argv)
