
import snakemake
import logging
import os
import re
import sys
import yaml
import argparse
from ete3 import NCBITaxa
from snakemake.utils import listfiles
from snakemake.workflow import expand


#needs to do:
#verification:
#parse the config file and determine that the fasta exist - should this be in the snakemake?
#determine that the buscos and input files are all where they should be
#make sure that taxdump is where it should be
#if no taxdump ete3 file, create it

#analysis:
#run just the alignments
#print just the alignment commands
#do the full filtering

#do everythign at once (let's start with this). this will be eukdetect run


def main(argv=sys.argv):

	parser = argparse.ArgumentParser(
		"eukdetect",
		description="Run full EukDetect pipeline."
		)

	parser.add_argument(
		"--configfile",
		type=str,
		action="store",
		dest="config",
		help="EukDetect config file (see default_configfile.yml"
		)

	args = parser.parse_args()

	#parse the yaml and make sure the files are all ok
	with open(args.config) as file:
		info = yaml.load(file, Loader=yaml.FullLoader)
		#here, if there is no readmin in the yaml file, write it
		#check read length is long enough

		if readlen < 70:
			#throw an error more gracefully than this
			sys.exit()

		if "readlenmin" not in info:
			readmin = str(round(info[readlen] ))
		print(info)


	#here we do the snakemake stuff




#def check_eukdb(db_path):


#def check_taxdump(db_path):

#def check_fastqs(fastq_path):