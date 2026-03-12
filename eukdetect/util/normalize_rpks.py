#!/usr/bin/env python3
"""
Create unified abundance table from multiple EukDetect samples

This script merges EukDetect eukfrac files into a single table with one row per species
and columns for each sample's RPKS, RPKSM (normalized), and EukFrac values.

Similar to MetaPhlAn's merge_metaphlan_tables.py output format.

Usage:
	python create_unified_table.py --eukfrac sample1_eukfrac.txt [sample2_eukfrac.txt ...] \\
	                               --library-sizes library_sizes.tsv \\
	                               --output unified_abundance.tsv

Library sizes file format (tab-separated):
	Column 1: sample name
	Column 2: total reads (integer)

Output format:
	TaxID  Name  Lineage  sample1_RPKS  sample1_RPKSM  sample1_EukFrac  sample2_RPKS  ...
"""

import argparse
import sys
import os
from pathlib import Path
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO
)
logger = logging.getLogger(__name__)


def parse_library_sizes(library_file):
	library_sizes = {}
	
	try:
		with open(library_file) as f:
			# Skip header
			f.readline()
			
			for line_num, line in enumerate(f, start=2):
				parts = line.strip().split('\t')
				
				if len(parts) < 2:
					logger.warning(f"Skipping line {line_num}: need at least 2 columns")
					continue
				
				sample = parts[0]
				
				try:
					total_reads = int(parts[1])
				except ValueError:
					raise ValueError(
						f"Error on line {line_num}: Column 2 must contain integers.\n"
						f"Sample: {sample}, Value: '{parts[1]}'"
					)
				
				if total_reads <= 0:
					logger.warning(f"Sample {sample} has invalid read count: {total_reads}")
					continue
				
				library_sizes[sample] = total_reads
		
		logger.info(f"Loaded library sizes for {len(library_sizes)} samples")
		return library_sizes
		
	except FileNotFoundError:
		logger.error(f"Library sizes file not found: {library_file}")
		sys.exit(1)
	except ValueError as e:
		logger.error(str(e))
		sys.exit(1)


def extract_sample_name(eukfrac_file):
	filename = Path(eukfrac_file).name
	
	suffixes_to_remove = [
		'_filtered_hits_eukfrac.txt',
		'_eukfrac.txt',
		'.txt'
	]
	
	sample_name = filename
	for suffix in suffixes_to_remove:
		if sample_name.endswith(suffix):
			sample_name = sample_name[:-len(suffix)]
			break
	
	return sample_name


def parse_eukfrac_file(eukfrac_file, library_sizes):
	sample_name = extract_sample_name(eukfrac_file)
	
	if sample_name not in library_sizes:
		logger.error(
			f"Sample '{sample_name}' not found in library sizes file.\n"
			f"Available samples: {', '.join(sorted(library_sizes.keys()))}"
		)
		sys.exit(1)
	
	total_reads = library_sizes[sample_name]
	#logger.info(f"Processing {sample_name}: {total_reads:,} total reads")
	
	taxon_data = {}
	
	try:
		with open(eukfrac_file) as f:
			# Read first line
			first_line = f.readline().strip()
			
			# Check for empty file indicators
			if not first_line:
				logger.warning(f"  Empty file - no data for {sample_name}")
				return sample_name, taxon_data
			
			if "Empty read count file" in first_line or "no aligned reads" in first_line.lower():
				logger.warning(f"  No aligned reads in {sample_name}")
				return sample_name, taxon_data
			
			if "No taxa passing filter" in first_line:
				logger.warning(f"  No taxa passed filters in {sample_name}")
				return sample_name, taxon_data
			
			# If we get here, first_line should be the header
			header = first_line.split('\t')
			
			# Check if this is actually a header (has expected columns)
			if 'Rank' not in header:
				# Might be a message or empty file
				logger.warning(f"  No valid data in {sample_name} (no header found)")
				return sample_name, taxon_data
			
			# Find required columns
			try:
				rank_idx = header.index('Rank')
				name_idx = header.index('Name')
				taxid_idx = header.index('TaxID')
				lineage_idx = header.index('Lineage')
				rpks_idx = header.index('RPKS')
				rel_abund_idx = header.index('Relative_abundance')
			except ValueError as e:
				logger.warning(f"  Missing required column in {sample_name}: {e}")
				return sample_name, taxon_data
			
			for line in f:
				# Check if line is a message
				if "Empty read count file" in line or "No taxa passing filter" in line:
					logger.warning(f"  No data in {sample_name}: {line.strip()}")
					return sample_name, taxon_data
				
				parts = line.strip().split('\t')
				
				if len(parts) != len(header):
					continue
				
				lineage = parts[lineage_idx]
				rank = parts[rank_idx]
				taxid = parts[taxid_idx]
				name = parts[name_idx]
				rpks_str = parts[rpks_idx]
				eukfrac_str = parts[rel_abund_idx]
				
				# Parse RPKS: NA for non-species ranks, numeric (including 0) for species
				if rank != 'species':
					rpks = 'NA'
					rpksm = 'NA'
				else:
					try:
						rpks = float(rpks_str) if rpks_str != 'NA' else 0.0
						rpksm = rpks / (total_reads / 1_000_000)
					except ValueError:
						rpks = 0.0
						rpksm = 0.0
				
				# Parse EukFrac
				try:
					eukfrac = float(eukfrac_str)
				except ValueError:
					eukfrac = 0.0
				
				taxon_data[lineage] = {
					'rank': rank,
					'name': name,
					'taxid': taxid,
					'rpks': rpks,
					'rpksm': rpksm,
					'eukfrac': eukfrac
				}
		
		if taxon_data:
			logger.info(f"  Found {len(taxon_data)} taxa")
		else:
			logger.warning(f"  No taxa found in {sample_name}")
		
		return sample_name, taxon_data
		
	except FileNotFoundError:
		logger.error(f"File not found: {eukfrac_file}")
		sys.exit(1)


def create_unified_table(eukfrac_files, library_sizes):

	#Collect data from all files
	sample_data = {}
	all_lineages = set()
	lineage_info = {}
	sample_order = []
	
	for eukfrac_file in eukfrac_files:
		if not os.path.exists(eukfrac_file):
			logger.warning(f"File not found: {eukfrac_file}, skipping")
			continue
		
		sample_name, taxon_data = parse_eukfrac_file(eukfrac_file, library_sizes)
		sample_order.append(sample_name)
		sample_data[sample_name] = taxon_data
		
		for lineage, data in taxon_data.items():
			all_lineages.add(lineage)
			if lineage not in lineage_info:
				lineage_info[lineage] = {
					'rank': data['rank'],
					'name': data['name'],
					'taxid': data['taxid']
				}
	
	logger.info(f"Total unique taxa across all samples: {len(all_lineages)}")
	
	return all_lineages, lineage_info, sample_data, sample_order


def write_unified_table(all_lineages, lineage_info, sample_data, sample_order, output_file):
	sorted_lineages = sorted(all_lineages)
	
	try:
		with open(output_file, 'w') as f:
			# Write header
			header_parts = ['Lineage', 'Rank', 'Name', 'TaxID']
			for sample in sample_order:
				header_parts.extend([
					f'{sample}_RPKS',
					f'{sample}_RPKSM',
					f'{sample}_EukFrac'
				])
			f.write('\t'.join(header_parts) + '\n')
			
			# Write data rows
			for lineage in sorted_lineages:
				info = lineage_info[lineage]
				row_parts = [lineage, info['rank'], info['name'], info['taxid']]
				
				for sample in sample_order:
					if sample in sample_data and lineage in sample_data[sample]:
						data = sample_data[sample][lineage]
						
						# Format RPKS and RPKSM (NA for non-species, numeric for species)
						if data['rpks'] == 'NA':
							rpks_str = 'NA'
							rpksm_str = 'NA'
						else:
							rpks_str = f"{data['rpks']:.6f}"
							rpksm_str = f"{data['rpksm']:.6f}"
						
						row_parts.extend([
							rpks_str,
							rpksm_str,
							f"{data['eukfrac']:.6f}"
						])
					else:
						# Taxon not found in this sample
						is_species = lineage_info[lineage]['rank'] == 'species'
						rpks_absent = '0.000000' if is_species else 'NA'
						row_parts.extend([rpks_absent, rpks_absent, '0.000000'])
				
				f.write('\t'.join(row_parts) + '\n')
		
		logger.info(f"Wrote unified table to {output_file}")
		logger.info(f"  {len(sorted_lineages)} taxa × {len(sample_order)} samples")
		
	except Exception as e:
		logger.error(f"Error writing output: {e}")
		sys.exit(1)


def main():
	parser = argparse.ArgumentParser(
		description="Create unified abundance table from EukDetect eukfrac files",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""
Examples:
  eukdetect-normalize \\
    --eukfrac results/*_filtered_hits_eukfrac.txt \\
    --library-sizes library_sizes.tsv \\
    --output unified_abundance.tsv

Output format:
  Lineage  Rank    Name      TaxID  sample1_RPKS  sample1_RPKSM  sample1_EukFrac  sample2_RPKS  ...
  phylum-Ascomycota  phylum  Ascomycota  4890  NA  NA  100.0  NA  ...
  phylum-Ascomycota|class-Saccharomycetes  class  Saccharomycetes  4891  NA  NA  100.0  NA  ...
  ...
  phylum-Ascomycota|...|species-Saccharomyces_cerevisiae  species  Saccharomyces cerevisiae  4932  0.0197  0.00394  100.0  ...

  - One row per taxon (all taxonomic levels: phylum, class, order, family, genus, species)
  - Three columns per sample: RPKS, RPKSM (normalized), EukFrac
  - RPKS and RPKSM are NA for non-species ranks
  - Zero values for taxa not detected in a sample
  - Taxonomic tree structure preserved (sorted by lineage string)
		"""
	)
	
	parser.add_argument(
		'--eukfrac',
		nargs='+',
		required=True,
		metavar='FILE',
		help='EukDetect eukfrac output files'
	)
	
	parser.add_argument(
		'--library-sizes',
		required=True,
		metavar='FILE',
		help='Tab-separated file: column 1 = sample names, column 2 = total reads'
	)
	
	parser.add_argument(
		'--output',
		required=True,
		metavar='FILE',
		help='Output file for unified table'
	)
	
	parser.add_argument(
		'--verbose',
		action='store_true',
		help='Enable verbose logging'
	)
	
	args = parser.parse_args()
	
	if args.verbose:
		logger.setLevel(logging.DEBUG)
	
	# Validate inputs
	if not args.eukfrac:
		logger.error("No eukfrac files provided")
		sys.exit(1)
	
	logger.info(f"Processing {len(args.eukfrac)} eukfrac file(s)")

	library_sizes = parse_library_sizes(args.library_sizes)
	
	if not library_sizes:
		logger.error("No library sizes loaded")
		sys.exit(1)

	all_lineages, lineage_info, sample_data, sample_order = create_unified_table(
		args.eukfrac, library_sizes
	)
	
	if not all_lineages:
		logger.error("No taxa found across all samples")
		sys.exit(1)
	

	write_unified_table(all_lineages, lineage_info, sample_data, sample_order, args.output)
	
	logger.info("-"*60)
	logger.info("Unified table creation complete!")
	logger.info(f"  Samples: {len(sample_order)}")
	logger.info(f"  Taxa: {len(all_lineages)}")
	logger.info(f"  Output: {args.output}")
	logger.info("-"*60)


if __name__ == "__main__":
	main()
