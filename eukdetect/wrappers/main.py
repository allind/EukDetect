from pathlib import Path

import argparse
import logging
import sys
import textwrap

from . import runall

logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO,
)

logger = logging.getLogger(__name__)

def create_parser():
	#main parser
	parser = argparse.ArgumentParser(
		prog = "eukdetect",
		description=textwrap.dedent("""\
			EukDetect: Detect and quantify abundance of eukaryotic microbes with metagenomic sequencing.

			Quick start:
				Single sample:
					eukdetect single -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o results/ -d database/

				Batch mode:
					eukdetect batch --samples samples.tsv -o results/ -d database/

		"""),

		formatter_class = argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument("--version", action='version', version='EukDetect v2.0.1')
	parser.add_argument("--verbose", "-v", action="store_true", help="verbose")

	#subparsers for commands
	subparsers = parser.add_subparsers(dest = "command", help = "Available commands", metavar="COMMAND")
	
	#Single sample command
	single_parser = subparsers.add_parser(
		"single", 
		help="Run single sample from command line", 
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	runall.parseargs_single(single_parser)

	#Batch command
	batch_parser = subparsers.add_parser(
		"batch",
		help="Run batch of samples from TSV file locally",
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	runall.parseargs_batch(batch_parser)
	
	return parser

def main():
	parser = create_parser()
	args = parser.parse_args()

	if args.verbose:
		#verbose
		logging.getLogger().setLevel(logging.DEBUG)
		logger.debug("Verbose mode")

	#print help if no args passed
	if not args.command:
		parser.print_help()
		sys.exit(0)

	#execute commands
	if args.command == "single":
		args.run_type = "single"  #Set for compatibility with execute()
		runall.execute(args)
	elif args.command == "batch":
		args.run_type = "batch"  #Set for compatibility with execute()
		runall.execute(args)
	else:
		parser.print_help()
		sys.exit(1)

if __name__ == "__main__":
	main()





