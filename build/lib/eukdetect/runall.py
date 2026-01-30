
#from ete3 import NCBITaxa

from Bio import SeqIO
from datetime import datetime
from pathlib import Path
import logging
import shlex
import gzip
import os
import sys
import yaml
import argparse
import subprocess
import textwrap


def main(argv=sys.argv):

	parser = argparse.ArgumentParser(
		description=textwrap.dedent("""\
			Run the EukDetect pipeline.

			Required arguments are --mode and --config.

			Modes:
				all: runs entire pipeline.

				aln: only runs alignment step of pipeline.

				filter: runs all downstream analysis of alignment. 
					Requires alignment files to already exist.

				printaln: prints alignment commands to alignment_commands.txt.


			For a schematic on EukDetect modes, see eukdetect_pipeline_schematic.pdf
			"""),
		formatter_class = argparse.RawDescriptionHelpFormatter
		)

	parser.add_argument(
		"--mode",
		type=str,
		action="store",
		dest="mode",
		help= "Run mode: choose from runall, aln, printaln, or filter.",
		required=True
		)

	parser.add_argument(
		"--configfile",
		type=str,
		action="store",
		dest="config",
		help="EukDetect config file (see example in default_configfile.yml)",
		required=True
		)

	parser.add_argument(
		"--cores",
		type=int,
		action="store",
		dest="cores",
		help="Number of cores to use for Snakemake analysis. Default is 1 core.",
		default=1
		)

	parser.add_argument(
		"--force",
		action="store_true",
		dest="force",
		help="Removes and rewrites existing analysis files in output directory.",
		default=False
		)

	log_level = logging.INFO
	logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=log_level,
    )
	options = parser.parse_args()


	logging.info("Parsing config file ...")

	#parse the yaml and make sure database files are where they should be

	with open(options.config) as file:
		config_info = yaml.load(file, Loader=yaml.FullLoader)
		#check read length is long enough

		if config_info["readlen"] < 75: #min readlength is 75 bp

			logging.error("Read length is below minimum 75 bp")
			exit(1)

		#check the database exists
		db_files = [config_info["database_prefix"] + ".1.bt2", 
					config_info["database_prefix"] + ".2.bt2", 
					config_info["database_prefix"] + ".3.bt2", 
					config_info["database_prefix"] + ".4.bt2", 
					config_info["database_prefix"] + ".rev.1.bt2", 
					config_info["database_prefix"] + ".rev.2.bt2",
					"taxa.sqlite",
					"taxa.sqlite.traverse.pkl",
					"specific_and_inherited_markers_per_taxid.txt",
					"busco_taxid_link.txt",
					"taxid_cumulativelength.txt"
					]

		if os.path.isdir(config_info["database_dir"]):
			seen_sql = False
			for f in os.listdir(config_info["database_dir"]):
				if f in db_files:
					db_files.remove(f)
				if f == "taxa.sqlite":
					seen_sql = True

			if len(db_files) > 0:
				logging.error("Error: missing files in database directory: " + 
					", ".join(db_files))
				exit(1)
		else:
			logging.error("Error: could not locate database directory.")
			exit(1)

	#check snakemake rules file exists
	snakefile = Path(config_info["eukdetect_dir"] + "/rules/eukdetect_eukfrac.rules")
	if not snakefile.exists():
		logging.error("Error: could not find /rules/eukdetect_eukfrac.rules in eukdetect_dir specified in configfile.")
		exit(1)

	#check required inputs and required outputs

	if options.mode == "runall":

		alncontain, alnmiss = check_aln_output(config_info)
		filtercontain, filtermiss = check_filter_output(config_info)
		fastacontain, fastamiss = check_fastas(config_info)
		fastalen = check_readlen(config_info)

		if len(fastamiss) > 0:
			logging.error("Error: missing input fasta files ...")
			for e in fastamiss:
				logging.error("Missing: " + e)
			exit(1)

		if len(fastalen) > 0:
			logging.error("Error: read length in these fastq files appears more than 10 basepairs different from config file read length.")
			logging.error("\n".join(fastalen))
			exit(1)

		if len(alncontain) > 0 or len(filtercontain) > 0:
			remain = alncontain + filtercontain
			if options.force:
				logging.info("Removing previous output files in output directory ...")
				for f in remain:
					os.remove(f)
			else:
				logging.error("Error: output files exist in output directory. Rerun EukDetect with --force to overwrite.")
				logging.error("\n".join(remain))
				exit(1)

		#check fastas exist

		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores)]

	elif options.mode == "aln":
		alncontain, alnmiss = check_aln_output(config_info)
		fastacontain, fastamiss = check_fastas(config_info)
		fastalen = check_readlen(config_info)

		if len(fastamiss) > 0:
			logging.error("Error: missing input fasta files ...")
			for e in fastamiss:
				logging.error("Missing: " + e)
			exit(1)
		if len(fastalen) > 0:
			logging.error("Error: read length in these fastq files appears more than 10 basepairs different from config file read length.")
			logging.error("\n".join(fastalen))
			exit(1)

		if len(alncontain) > 0:
			if options.force:
				logging.info("Removing previous output files in output directory ...")
			else:
				logging.error("Error: output files exist in output directory. Rerun EukDetect with --force to overwrite.")
				logging.error("\n".join(alncontain))
				exit(1)

		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores),
						  "aln"
						  ]

	elif options.mode == "printaln":
		fastacontain, fastamiss = check_fastas(config_info)
		fastalen = check_readlen(config_info)

		if len(fastamiss) > 0:
			logging.error("Error: missing input fasta files ...")
			for e in fastamiss:
				logging.error("Missing: " + e)
			exit(1)
		if len(fastalen) > 0:
			logging.error("Error: read length in these fastq files appears more than 10 basepairs different from config file read length.")
			logging.error("\n".join(fastalen))
			exit(1)

		#remove old alignment commands file if it exists
		if os.path.isfile(config_info["output_dir"] + "/alignment_commands.txt"):
			if options.force:
				logging.info("Removing previous alignment command file in output directory ...")
				os.remove(config_info["output_dir"] + "/alignment_commands.txt")
			else:
				logging.error("Output directory already has alignment_commands.txt file. Rerun EukDetect with --force to overwrite.")
				exit(1)



		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores),
						  "printaln"
						  ]



	elif options.mode == "filter":
		#make sure alignment files are where they are supposed to be
		alncontain, alnmiss = check_aln_output(config_info)

		if len(alnmiss) > 0:
			logging.error("Missing required files from input step:")
			for e in alnmiss:
				logging.error("Missing: " + e)
			exit(1)


		#check for existing filter files
		filtercontain, filtermiss = check_filter_output(config_info)
		if len(filtercontain) > 0:
			if options.force:
				logging.info("Removing previous output files in output directory ...")
				for e in filtercontain:
					os.remove(e)
			else:
				logging.error("Error: output files exist in output directory. Rerun EukDetect with --force to overwrite.")
				for e in filtercontain:
					logging.error(e)
				exit(1)


		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores),
						  "filter"
						  ]
	else:
		logging.error("Error: incorrect mode specified.")
		exit(1)


	#call snakemake

	logging.info("Running: " + " ".join(snakemake_args))
	timestamp = datetime.timestamp(datetime.now())
	logging.info("Redirecting snakemake output to snakemake_" + str(timestamp) + ".log")

	snakelog = open("snakemake_" + str(timestamp) + ".log", "w")
	cmd = subprocess.run(snakemake_args, stdout=snakelog, stderr=snakelog)
	logging.info("Snakemake finished")

	if cmd.returncode == 0:
		#check correct output exists

		if options.mode == "printaln":
			#if printaln replace the brackets with single quotes that snakemake improperly evalutes. would like to find alt solution	
			if os.path.isfile(config_info["output_dir"] + "/alignment_commands.txt"):
				newcmds = []
				for line in open(config_info["output_dir"] + "/alignment_commands.txt"):
					newcmds.append(line.replace("{}", "'"))
				os.remove(config_info["output_dir"] + "/alignment_commands.txt")
				dest = open(config_info["output_dir"] + "/alignment_commands.txt", 'w')
				for c in newcmds:
					dest.write(c)
				dest.close()
				logging.info("Snakemake pipeline created all files.")
				exit(0)

			else:
				logging.error("Something went wrong with printing alignment commands. Check logs.")
				exit(1)

		elif options.mode == "aln":
			contains, missing = check_aln_output(config_info)
			if len(missing) > 0:
				logging.error("Something went wrong with alignment. Files are missing. Check logs.")
				for e in missing:
					logging.error("Missing: " + e)
				exit(1)
			else:
				logging.info("Snakemake pipeline created all files.")
				exit(0)

		elif options.mode == "runall":
			acontains, amissing = check_aln_output(config_info)
			fcontains, fmissing = check_filter_output(config_info)
			contains = acontains + fcontains
			missing = amissing + fmissing
			if len(missing) > 0:
				logging.error("Something went wrong. Files are missing. Check logs.")
				for e in missing:
					logging.error("Missing: " + e)
				exit(1)
			else:
				logging.info("Snakemake pipeline created all files.")
				exit(0)




		elif options.mode == "filter":
			contains, missing = check_filter_output(config_info)
			if len(missing) > 0:
				logging.error("Something went wrong with alignment. Files are missing. Check snakemake logs.")
				for e in missing:
					logging.error("Missing: " + e)
				exit(1)
			else:
				logging.info("Snakemake pipeline created all files.")
				exit(0)



	if cmd.returncode != 0:
		logging.error("Something went wrong with snakemake. Check log files")
		exit(cmd.returncode)




#output checking functions

def check_aln_output(config_info):

	samples = config_info["samples"].keys()
	exists = []
	outdir = config_info["output_dir"]
	contains = []
	missing = []

	for sample in samples:
		if os.path.isfile(outdir + '/aln/' + sample + '_aln_q30_lenfilter.sorted.bam'):
			contains.append(outdir + '/aln/' + sample + '_aln_q30_lenfilter.sorted.bam')
		else:
			missing.append(outdir + '/aln/' + sample + '_aln_q30_lenfilter.sorted.bam')
	return contains, missing

def check_filter_output(config_info):

	samples = config_info["samples"].keys()
	exists = []
	outdir = config_info["output_dir"]
	missing = []
	contains = []

	for sample in samples:
		f1 = outdir + '/filtering/' + sample + "_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam"
		f2 = outdir + '/filtering/' + sample + '_aln_q30_lenfilter_complexityfilter_dupfilter.sorted.bam.bai'
		f3 = outdir + '/filtering/' + sample + "_read_counts_and_mismatches.txt"
		f4 = outdir + '/filtering/' + sample + "_all_hits_table.txt"
		f5 = outdir + '/' + sample + "_filtered_hits_table.txt"
		f6 = outdir + '/' + sample + "_filtered_hits_eukfrac.txt"

		filelist = [f1, f2, f3, f4, f5, f6]

		for f in filelist:
			if os.path.isfile(f):
				contains.append(f)
			else:
				missing.append(f)
	return contains, missing

def check_fastas(config_info):
	samples = config_info["samples"].keys()
	contains = []
	missing = []
	path = config_info["fq_dir"] + "/"
	if config_info["paired_end"]:
		for sample in samples:
			fwd = path + sample + config_info["fwd_suffix"]
			rev = path + sample + config_info["rev_suffix"]
			if os.path.isfile(fwd):
				contains.append(fwd)
			else:
				missing.append(fwd)

			if os.path.isfile(rev):
				contains.append(rev)
			else:
				missing.append(rev)
	else:
		for sample in samples:
			read = path + sample + config_info["se_suffix"]
			if os.path.isfile(read):
				contains.append(read)
			else:
				contains.append(read)

	return contains, missing

def check_readlen(config_info):
	readlen = config_info["readlen"]
	samples = config_info["samples"].keys()
	path = config_info["fq_dir"] + "/"

	bad = []
	if config_info["paired_end"]:
		for sample in samples:
			fwd = path + sample + config_info["fwd_suffix"]
			rev = path + sample + config_info["rev_suffix"]
		if os.path.isfile(fwd):
			#check readlen
			counter = 0
			bases = 0
			if ".gz" in fwd:
				with gzip.open(fwd, "rt") as handle:
					for record in SeqIO.parse(handle, "fastq"):
						if counter > 10000:
							break
						counter += 1
						bases += len(str(record.seq))
			else:
				for record in SeqIO.parse(fwd, "fastq"):
					if counter > 10000:
						break
					counter += 1
					bases += len(str(record.seq))

			file_readlen = int(bases / counter)
			if abs(file_readlen - readlen) > 10:
				bad.append(fwd)

		if os.path.isfile(rev):
			#check readlen
			counter = 0
			bases = 0
			if ".gz" in rev:
				
				with gzip.open(rev, "rt") as handle:
					for record in SeqIO.parse(handle, "fastq"):
						if counter > 10000:
							break
						counter += 1
						bases += len(str(record.seq))
			else:
				for record in SeqIO.parse(rev, "fastq"):
					if counter > 10000:
						break
					counter += 1
					bases += len(str(record.seq))

			file_readlen = int(bases / counter)
			if abs(file_readlen - readlen) > 10:
				bad.append(rev)

	else:
		for sample in samples:
			read = path + sample + config_info["se_suffix"]
			if os.path.isfile(read):
				counter = 0
				bases = 0
				if ".gz" in read:
					with gzip.open(read, "rt") as handle:
						for record in SeqIO.parse(handle, "fastq"):
							if counter > 10000:
								break
							counter += 1
							bases += len(str(record.seq))
				else:
					for record in SeqIO.parse(read, "fastq"):
						if counter > 10000:
							break
						counter += 1
						bases += len(str(record.seq))

				file_readlen = int(bases / counter)
				if abs(file_readlen - readlen) > 10:
					bad.append(read)
	return bad

