
#from ete3 import NCBITaxa

from datetime import datetime
from pathlib import Path
#from check_output import check_aln_output, check_filter_output
#from eukdetect import check_output as chk
import logging
import os
import sys
import yaml
import argparse
import subprocess


def main(argv=sys.argv):

	parser = argparse.ArgumentParser(description="Run full EukDetect pipeline.")

	parser.add_argument(
		"--mode",
		type=str,
		action="store",
		dest="mode",
		help="Run mode: choose from runall, aln, alncmd, or filter.",
		required=True
		)

	parser.add_argument(
		"--configfile",
		type=str,
		action="store",
		dest="config",
		help="EukDetect config file (see default_configfile.yml)",
		required=True
		)

	parser.add_argument(
		"--cores",
		type=int,
		action="store",
		dest="cores",
		help="Number of cores to use",
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

		#check the database
		db_files = ["all_buscos_v4.fna.1.bt2", 
					"all_buscos_v4.fna.2.bt2", 
					"all_buscos_v4.fna.3.bt2", 
					"all_buscos_v4.fna.4.bt2", 
					"all_buscos_v4.fna.rev.1.bt2", 
					"all_buscos_v4.fna.rev.2.bt2",
					"taxa.sqlite"]

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

	#here we do the snakemake stuff
	snakefile = Path(config_info["eukdetect_dir"] + "/rules/eukdetect.rules")
	if not snakefile.exists():
		logging.error("Error: could not find /rules/runall.rules in eukdetect_dir specified in configfile.")
		exit(1)

	if options.mode == "runall":
		alnmiss, alncontain = check_aln_output(config_info)
		filtermiss, filtercontain = check_filter_output(config_info)

		if len(alncontain) > 0 or len(filtercontain) > 0:
			remain = alncontain + filtercontain
			if options.force:
				logging.info("Removing previous output files in output directory ...")
				for f in remain:
					os.remove(f)
			else:
				logging.error("Error: output files exist in output directory...")
				logging.error("\n".join(remain))
				exit(1)

		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores)]

	elif options.mode == "aln":
		alnmiss, alncontain = check_aln_output(config_info)
		if len(alncontain) > 0:
			if options.force:
				logging.info("Removing previous output files in output directory ...")
			else:
				logging.error("Error: output files exist in output directory...")
				logging.error("\n".join(alncontain))
				exit(1)

		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores),
						  "runaln"
						  ]

	elif options.mode == "alncmd":

		#remove old alignment commands file if it exists
		if os.path.isfile(config_info["output_dir"] + "/alignment_commands.txt"):
			if options.force:
				logging.info("Removing previous alignment command file in output directory ...")
				os.remove(config_info["output_dir"] + "/alignment_commands.txt")
			else:
				logging.error("Output directory already has alignment_commands.txt file.")
				exit(1)

		snakemake_args = ['snakemake', 
						  '--snakefile', 
						  str(snakefile), 
						  '--configfile',
						  options.config,
						  '--cores',
						  str(options.cores),
						  "alncall"
						  ]



	elif options.mode == "filter":
		#make sure alignment files are where they are supposed to be
		filtercontain, filtermiss = check_filter_output(config_info)
		if len(filtercontain) > 0:
			if options.force:
				logging.info("Removing previous output files in output directory ...")
				for e in filtercontain:
					os.remove(e)
			else:
				logging.error("Error: output files exist in output directory. Shutting down.")
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
	logging.info("Snakemake complete")


	
	
	#if alncmd replace the brackets with single quotes that snakemake improperly evalutes. would like to find alt solution

	if cmd.returncode == 0:

		if options.mode == "alncmd":	
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
				logging.error("Something went wrong with printing alignment commands. Check logs")
				exit(1)

		elif options.mode == "aln":
			contains, missing = check_aln_output(config_info)
			if len(missing) > 0:
				logging.error("Something went wrong with alignment. Files are missing. Check snakemake logs.")
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
				logging.error("Something went wrong. Files are missing. Check snakemake logs.")
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



	if cmd.returncode != 1:
		logging.error("Something went wrong with snakemake. Check log files")
		exit(cmd.returncode)




#output checking function

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
		f4 = outdir + '/filtering/' + sample + "_hit_taxonomy.txt"
		f5 = outdir + '/filtering/' + sample + "_all_stats_per_taxid.txt"
		f6 = outdir + '/' + sample + "_hit_taxonomy_filterpass.txt"
		f7 = outdir + '/' + sample + "_stats_per_filtered_taxid.txt"

		filelist = [f1, f2, f3, f4, f5, f6, f7]

		for f in filelist:
			if os.path.isfile(f):
				contains.append(f)
			else:
				missing.append(f)
	return contains, missing
