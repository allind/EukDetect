
from ete3 import NCBITaxa
from pathlib import Path
import logging
import os
import sys
import yaml
import argparse
import subprocess



#needs to do:
#verification:
#parse the config file and determine that the fasta exist - should this be in the snakemake?
#determine that the buscos and input files are all where they should be
#make sure that taxdump is where it should be
#if no taxdump ete3 file, create it. actually, I can't really do this. This might be an issue.

#analysis:
#run just the alignments
#print just the alignment commands
#do the full filtering

#do everythign at once (let's start with this). this will be eukdetect run


def main(argv=sys.argv):

	parser = argparse.ArgumentParser(description="Run full EukDetect pipeline.")

	parser.add_argument(
		"--configfile",
		type=str,
		action="store",
		dest="config",
		help="EukDetect config file (see default_configfile.yml",
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

	log_level = logging.INFO
	logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=log_level,
    )

	logging.info("Parsing config file ...")

	options = parser.parse_args()

	#parse the yaml and make sure database files are ok

	with open(options.config) as file:
		config_info = yaml.load(file, Loader=yaml.FullLoader)
		#check read length is long enough

		if config_info["readlen"] < 75:
			#throw an error more gracefully than this
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
			logging.error("Error: ould not locate database directory.")
			exit(1)
		
	#here we do the snakemake stuff
	snakefile = Path(config_info["eukdetect_dir"] + "/rules/runall.rules")
	if not snakefile.exists():
		logging.error("Error: could not find /rules/runall.rules in eukdetect_dir specified in configfile.")
		exit(1)

	snakemake_args = ['snakemake', 
					  '--snakefile', 
					  str(snakefile), 
					  '--configfile',
					  options.config,
					  '--cores',
					  str(options.cores)]

	logging.info("Running: " + " ".join(snakemake_args))
	cmd = subprocess.run(snakemake_args)
	exit(cmd.returncode)
