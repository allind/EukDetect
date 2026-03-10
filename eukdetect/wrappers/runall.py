#from ete3 import NCBITaxa
#from Bio import SeqIO
#from datetime import datetime
from pathlib import Path
from typing import Optional, List, Tuple

import logging
import shlex
#import gzip
import os
import sys
import yaml
import argparse
#import subprocess
import textwrap
import csv
import re

from ..util.execute import SnakemakeExecutor
from ..util.build_config import ConfigBuilder
from ..util.validate import validate_inputs

logger = logging.getLogger(__name__)


def validate_cores(value):
	#Check that cores seems reasonable
	try:
		ivalue = int(value)
	except ValueError:
		raise argparse.ArgumentTypeError(f"--cores must be an integer, got '{value}'")
	
	if ivalue < 1:
		raise argparse.ArgumentTypeError("--cores must be at least 1")
	
	max_cores = os.cpu_count()
	if max_cores and ivalue > max_cores * 2:
		logger.warning(
			f"Requested {ivalue} cores but system has {max_cores} physical cores. "
			f"This may cause performance issues."
		)
	
	return ivalue


def validate_sample_name(name):
	#Check for valid characters
	if not re.match(r'^[a-zA-Z0-9_.-]+$', name):
		raise ValueError(
			f"Invalid sample name: '{name}'\n"
			f"Sample names can only contain letters, numbers, underscore, dash, period"
		)
	
	#Check it doesn't start with problematic characters
	if name.startswith('.') or name.startswith('-'):
		raise ValueError(f"Sample name cannot start with '.' or '-': '{name}'")
	
	#Check for path traversal attempts
	if '..' in name or '/' in name or '\\' in name:
		raise ValueError(f"Sample name cannot contain path separators: '{name}'")
	
	return name



def parseargs_single(parser):
	
	parser.description = textwrap.dedent("""\
		Run EukDetect on a single sample from command line.
		
		Modes:
			all       - Run full pipeline (alignment, filtering, analysis). Default.
			aln       - Run alignment step only.
			analyze   - Run filtering and analysis only (requires existing alignment).
			printaln  - Generate alignment commands file only.
		
		Example:
			eukdetect single -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \\
				-n mysample -o results/ -d /path/to/database --cores 8
	""")
	
	inputs = parser.add_argument_group("Inputs")
	
	inputs.add_argument(
		"-1", "--reads1",
		dest="reads1",
		required=True,
		metavar="PATH",
		help="Forward reads file for paired-end or single-end file (required)"
	)
	
	inputs.add_argument(
		"-2", "--reads2",
		dest="reads2",
		metavar="PATH",
		help="Reverse reads file for paired-end sequencing"
	)
	
	inputs.add_argument(
		"-n", "--name",
		dest="sample_name",
		metavar="NAME",
		help="Sample name (auto-detect from filename if not provided)"
	)
	
	outputs = parser.add_argument_group("Outputs")
	
	outputs.add_argument(
		"--outdir", "-o",
		dest="output",
		required=True,
		metavar="PATH",
		help="Output directory for results (required)"
	)
	
	dbopt = parser.add_argument_group("Database options")
	
	dbopt.add_argument(
		"--database", "-d",
		required=True,
		metavar="PATH",
		help="Path to eukdetect database directory (required)"
	)
	
	dbopt.add_argument(
		"--database-prefix",
		dest="database_prefix",
		default="eukdb",
		metavar="NAME",
		help="Database prefix name (default: eukdb)"
	)
	
	execs = parser.add_argument_group("Execution")
	
	execs.add_argument(
		"--mode", "-m",
		choices=['all', 'aln', 'analyze', 'printaln'],
		default='all',
		help="Analysis mode (default: all)"
	)
	
	execs.add_argument(
		"--cores", "-c",
		type=validate_cores,
		default=1,
		metavar="N",
		help="Number of CPU threads for bowtie2 (default: 1)"
	)
	
	execs.add_argument(
		"--readlen",
		type=int,
		metavar="N",
		help="Read length (auto-detect from first file if not provided)"
	)
	
	execs.add_argument(
		"--force",
		action="store_true",
		help="Overwrite existing files"
	)
	
	execs.add_argument(
		"--dry-run",
		action="store_true",
		help="Print commands without executing"
	)


def parseargs_batch(parser):	
	parser.description = textwrap.dedent("""\
		Run EukDetect on multiple samples locally.
		
		For execution on a compute cluster, use 'eukdetect single' paired with a job array.
		
		Analysis modes:
			all       - Run full pipeline (alignment, filtering, analysis). Default.
			aln       - Run alignment step only.
			analyze   - Run filtering and analysis only (requires existing alignment).
			printaln  - Generate alignment commands file only.
		
		Examples:
			#Run 10 samples in parallel locally
			eukdetect batch --samples samples.tsv -o results/ -d /path/to/database --cores 10
			
			#Generate alignment commands (for manual execution or cluster submission)
			eukdetect batch --samples samples.tsv -o results/ -d /path/to/database --mode printaln
		
		samples.tsv format (tab-separated):
			sample_name    reads1_path    reads2_path
			sample1        /path/to/sample1_R1.fq.gz    /path/to/sample1_R2.fq.gz
			sample2        /path/to/sample2_R1.fq.gz    /path/to/sample2_R2.fq.gz
	""")
	
	inputs = parser.add_argument_group("Inputs")
	
	inputs.add_argument(
		"--samples", "-s",
		required=True,
		metavar="PATH",
		help="Tab-separated file with sample information (required). Columns: sample_name, reads1, reads2"
	)
	
	inputs.add_argument(
		"--config",
		metavar="PATH",
		help="Config file in YAML format (optional, can override with other flags)"
	)
	
	outputs = parser.add_argument_group("Outputs")
	
	outputs.add_argument(
		"--outdir", "-o",
		dest="output",
		required=True,
		metavar="PATH",
		help="Output directory for results (required)"
	)
	
	outputs.add_argument(
		"--configfile-out",
		metavar="PATH",
		help="Write generated config file to this path (default: auto-generated name)"
	)
	
	dbopt = parser.add_argument_group("Database options")
	
	dbopt.add_argument(
		"--database", "-d",
		required=True,
		metavar="PATH",
		help="Path to eukdetect database directory (required)"
	)
	
	dbopt.add_argument(
		"--database-prefix",
		dest="database_prefix",
		default="eukdb",
		metavar="NAME",
		help="Database prefix name (default: eukdb)"
	)
	
	execs = parser.add_argument_group("Execution")
	
	execs.add_argument(
		"--mode", "-m",
		choices=['all', 'aln', 'analyze', 'printaln'],
		default='all',
		help="Analysis mode (default: all)"
	)
	
	execs.add_argument(
		"--cores", "-c",
		type=validate_cores,
		default=1,
		metavar="N",
		help="Number of samples to run in parallel on this machine (default: 1)."
	)
	
	execs.add_argument(
		"--readlen",
		type=int,
		metavar="N",
		help="Read length (auto-detect from reads1 file if not provided)"
	)
	
	execs.add_argument(
		"--force",
		action="store_true",
		help="Overwrite existing files"
	)
	
	execs.add_argument(
		"--dry-run",
		action="store_true",
		help="Show commands that will be run without running anything."
	)


def execute(args):
	try:
		database = args.database
		
		#Set bowtie2 cores and snakemake cores
		if args.run_type == 'single':
			bowtie2_cores = args.cores  #All cores go to bowtie2 threading
			snakemake_cores = args.cores
			logger.info(f"Running in SINGLE mode: bowtie2 will use {bowtie2_cores} threads")
		else:  #batch - local execution only
			bowtie2_cores = 1  #One thread per sample
			snakemake_cores = args.cores  #Number of parallel samples
			logger.info(f"Running in BATCH mode (local execution):")
			logger.info(f"  - {snakemake_cores} samples in parallel on this machine")
			logger.info(f"  - 1 thread per bowtie2 instance")
			logger.info(f"  - Total cores used: {snakemake_cores}")

		if hasattr(args, 'config') and args.config:
			logger.info(f"Using config file: {args.config}")
			config_dict = _load_config_file(args.config)

			#command line options can override config
			if database:
				config_dict["database_dir"] = database
			if args.output:
				config_dict["output_dir"] = args.output
			
			#Update bowtie2_cores in config
			config_dict["bowtie2_cores"] = bowtie2_cores
		else:
			logger.info("Building config from command-line arguments")
			
			if args.run_type == 'single':
				if not args.reads1:
					raise ValueError("Single mode requires -1/--reads1 argument")
				
				reads1 = args.reads1
				reads2 = args.reads2 if hasattr(args, 'reads2') else None
				sample_name = args.sample_name if hasattr(args, 'sample_name') and args.sample_name else None
				
				sample_dict = _parse_samples(
					reads1=[reads1],  #Convert to list
					reads2=[reads2] if reads2 else [],
					sample_name=[sample_name] if sample_name else [],
					samples_file=None
				)
			else:

				if not args.samples:
					raise ValueError("Batch mode requires --samples argument")
				
				#Check that samples file exists
				if not Path(args.samples).exists():
					raise FileNotFoundError(
						f"Samples file not found: {args.samples}\n"
						f"Please provide a valid tab-separated file with columns: sample_name, reads1, reads2"
					)
				
				sample_dict = _parse_samples(
					reads1=[],
					reads2=[],
					sample_name=[],
					samples_file=args.samples
				)

			if not database:
				raise ValueError(
					"Database directory required. Use --database."
					)
			if not sample_dict:
				raise ValueError(
					"No samples specified."
					)

			#determine if paired end or single end reads
			paired_end = any("reads2" in info for info in sample_dict.values())

			#create config dict
			config_builder = ConfigBuilder(
				samples=sample_dict,
				output_dir=args.output,
				database_dir=database,
				database_prefix=args.database_prefix,
				paired_end=paired_end,
				readlen=args.readlen,
				bowtie2_cores=bowtie2_cores,  #NEW
			)

			config_dict = config_builder.build()
			logger.info("Validating input")
			validate_inputs(config_dict, mode=args.mode, force=args.force)

			if hasattr(args, 'configfile_out') and args.configfile_out:
				config_path = Path(args.configfile_out)
			else:

				output_dir = Path(args.output)
				config_dir = output_dir / "configs"
				config_dir.mkdir(parents=True, exist_ok=True)
				
				if args.run_type == 'single':
					# Single mode: use sample name + run number
					sample_names = list(config_dict['samples'].keys())
					if sample_names:
						sample_name = sample_names[0]
						run_num = 1
						while (config_dir / f"config_{sample_name}_{run_num}.yml").exists():
							run_num += 1
						config_path = config_dir / f"config_{sample_name}_{run_num}.yml"
					else:
						#Fallback if no sample name found
						config_path = config_dir / "config_single.yml"
				else:
					#Batch mode: use run number
					run_num = 1
					while (config_dir / f"config_run{run_num}.yml").exists():
						run_num += 1
					config_path = config_dir / f"config_run{run_num}.yml"
			

			logger.info(f"Writing config to {config_path}")
			Path(config_path).parent.mkdir(parents=True, exist_ok=True)
			try:
				with open(config_path, "w") as f:
					yaml.dump(config_dict, f, default_flow_style=False)
			except (IOError, OSError) as e:
				raise IOError(
					f"Failed to write config file to {config_path}: {e}\n"
					f"Please check disk space and permissions."
				)

		#Warn if alignment files already exist when running in all or aln mode
		if args.mode in ('all', 'aln'):
			_warn_existing_alignments(config_dict)

		#Execute snakemake
		logger.info(f"Running EukDetect in {args.mode} mode")
		
		#Print dry-run info if applicable
		is_dry_run = hasattr(args, 'dry_run') and args.dry_run
		if is_dry_run and args.run_type == 'batch':
			num_samples = len(config_dict.get('samples', {}))
			logger.info("\n" + "*"*70)
			logger.info("DRY RUN - Nothing will be executed")
			logger.info("*"*70)
			logger.info(f"Number of samples: {num_samples}")
			logger.info(f"Samples that would run in parallel: {snakemake_cores}")
			logger.info(f"Mode: {args.mode}")
			logger.info("*"*70 + "\n")

		executor = SnakemakeExecutor(
			config_dict=config_dict,
			mode=args.mode,
			cores=snakemake_cores,
			force=args.force,
			dry_run=is_dry_run,
			)

		success = executor.run()

		if success:
			if is_dry_run:
				logger.info("\n" + "*"*70)
				logger.info("DRY RUN COMPLETED")
				logger.info("*"*70)
				if args.run_type == 'batch':
					logger.info("Dry run complete. To execute, re-run and remove --dry-run:")
					logger.info(f"  eukdetect batch --samples {args.samples} -o {args.output} -d {database} --cores {args.cores}")
				logger.info("*"*70 + "\n")
			else:
				logger.info("EukDetect completed.")
				_print_output_summary(config_dict, args.mode)
			sys.exit(0)
		else:
			logger.error("Eukdetect failed. Check logs for details.")
			sys.exit(1)

	except Exception as e:
		logger.error(f"Error: {e}")
		if hasattr(args, 'verbose') and args.verbose:
			raise
		sys.exit(1)


def _load_config_file(config_path: str) -> dict:
	with open(config_path) as f:
		config = yaml.load(f, Loader=yaml.FullLoader)
	return config

def _parse_samples(
	reads1: List[str],
	reads2: List[str],
	sample_name: List[str],
	samples_file: Optional[str],
) -> dict:

	#returns {sample_name: {'reads1': path, 'reads2': path}}. Reads2 omitted if they are not present and sample is assumed to be single end

	samples = {}
	if samples_file:
		#parsing from samples tsv
		logger.info(f"Loading samples from {samples_file}")
		with open(samples_file) as f:
			reader = csv.DictReader(f, delimiter="\t")
			for row in reader:
				if "sample_name" not in row or "reads1" not in row:
					raise ValueError(
						"Sample file must have columns: sample_name, reads1, [reads2]")
				sample = row["sample_name"]
				
				validate_sample_name(sample)
				
				samples[sample] = {"reads1": row["reads1"]}
				if "reads2" in row and row["reads2"]:
					samples[sample]["reads2"] = row["reads2"]
	elif reads1:
		logger.info(f"Parsing sample from command line input.")
		#Process reads1 file
		for idx, r1 in enumerate(reads1):
			if sample_name and idx < len(sample_name) and sample_name[idx]:
				#Use provided sample name
				name = sample_name[idx]
			else:
				#Auto-detect from filename - this is hacky
				name = Path(r1).name  

				extensions_to_remove = [
					".fastq.gz", ".fq.gz", ".fastq", ".fq",
					".FASTQ.GZ", ".FQ.GZ", ".FASTQ", ".FQ"
				]
				for ext in extensions_to_remove:
					if name.endswith(ext):
						name = name[:-len(ext)]
						break  #Only remove first matching extension
				
				#Remove read pair indicators
				pair_suffixes = ["_R1", "_R2", "_1", "_2", "_fwd", "_rev", "_F", "_R", "_FWD", "_REV"]
				for suffix in pair_suffixes:
					if name.endswith(suffix):
						name = name[:-len(suffix)]
						break  #Only remove first matching suffix
			
			#Validate sample name
			validate_sample_name(name)
			
			samples[name] = {"reads1": r1}
			if reads2 and idx < len(reads2) and reads2[idx]:
				samples[name]["reads2"] = reads2[idx]

	return samples

def _warn_existing_alignments(config_dict: dict):
	output_dir = Path(config_dict["output_dir"])
	samples = config_dict["samples"].keys()

	existing = [
		s for s in samples
		if (output_dir / "aln" / f"{s}_aln_q10_lenfilter.sorted.bam").exists()
	]

	if existing:
		logger.warning(
			f"\n{'!'*60}\n"
			f"WARNING: Alignment files already exist for {len(existing)}/{len(list(samples))} sample(s):\n"
			+ "\n".join(f"  - {s}" for s in existing)
			+ f"\n\nExisting alignments will NOT be re-run (Snakemake treats them as up-to-date)."
			f"\nTo redo alignments, delete the BAM files first."
			f"\nTo skip alignment entirely, use --mode analyze."
			f"\n{'!'*60}\n"
		)


def _print_output_summary(config_dict: dict, mode: str):
	output_dir = Path(config_dict["output_dir"])
	samples = config_dict["samples"].keys()
	logger.info("\n" + "*"*60)
	logger.info("Output summary")
	logger.info("*"*60)
	if mode in ["all", "analyze"]:
		logger.info(f"\nResults directory: {output_dir}")
		logger.info("\nFiltered output files:")
		for sample in samples:
			hit_table = output_dir / f"{sample}_filtered_hits_table.txt"
			eukfrac = output_dir / f"{sample}_filtered_hits_eukfrac.txt"
			logger.info(f"  • {sample}:")
			logger.info(f"	  - {hit_table}")
			logger.info(f"	  - {eukfrac}")

	elif mode == "aln":
		logger.info(f"\nAlignment output directory: {output_dir / 'aln'}")
		logger.info("\nAlignment bam files created for:")
		for sample in samples:
			logger.info(f"  • {sample}")

	elif mode == "printaln":
		cmd_file = output_dir / "alignment_commands.txt"
		logger.info(f"\nAlignment commands written to: {cmd_file}")
	logger.info("\n" + "*"*60)
