from Bio import SeqIO
from pathlib import Path
from typing import List, Tuple

import logging
import os
import gzip

logger = logging.getLogger(__name__)



def validate_inputs(config: dict, mode: str, force: bool) -> None:

	#Check output directory
	output_dir = Path(config["output_dir"])
	output_dir.mkdir(parents=True, exist_ok=True)
	
	#Check database
	check_database(config)
	
	#Check fastq files exist
	if mode in ["all", "aln", "printaln"]:
		check_fastq_files(config)
		check_readlengths(config)
	
	#Check for existing output files
	if mode in ["all", "aln"]:
		check_alignment_outputs(config, force)
	
	if mode in ["all", "analyze"]:
		check_filter_outputs(config, force)
	
	if mode == "analyze":
		check_alignment_inputs(config)
	
	logger.info("Validations passed")


def check_database(config: dict) -> None:
	db_dir = Path(config["database_dir"])
	db_prefix = config["database_prefix"]
	
	if not db_dir.is_dir():
		raise ValueError(f"Database directory not found: {db_dir}")
	
	# Check non-bowtie2 files
	required_files = [
		"taxa.sqlite",
		"taxa.sqlite.traverse.pkl",
		"busco_taxid_genome_link.txt",
		"taxid_and_genome_cumulativelength.txt",
	]
	
	missing = []
	for filename in required_files:
		if not (db_dir / filename).exists():
			missing.append(filename)
	
	# Check for bowtie2 index - either .bt2l (large) OR .bt2 (standard)
	bt2_variants = ['.bt2l', '.bt2']
	bt2_suffixes = ['.1', '.2', '.3', '.4', '.rev.1', '.rev.2']
	
	found_variant = None
	for variant in bt2_variants:
		if (db_dir / f"{db_prefix}.1{variant}").exists():
			found_variant = variant
			break
	
	if not found_variant:
		raise ValueError(
			f"Bowtie2 index not found for prefix '{db_prefix}' in {db_dir}\n"
			f"Expected either {db_prefix}.1.bt2l or {db_prefix}.1.bt2"
		)
	
	logger.debug(f"Found bowtie2 index with extension: {found_variant}")
	
	# Check all required bt2 files with the detected variant
	for suffix in bt2_suffixes:
		bt2_file = f"{db_prefix}{suffix}{found_variant}"
		if not (db_dir / bt2_file).exists():
			missing.append(bt2_file)
	
	if missing:
		raise ValueError(
			f"Missing database files in {db_dir}:\n  " + "\n  ".join(missing)
		)


	ani_path = config.get("ani_file") or str(db_dir / "ani.tsv")
	if not Path(ani_path).exists():
		raise ValueError(
			f"ANI file not found: {ani_path}\n"
			f"Expected three tab-separated columns with no header: "
			f"genome_1, genome_2, ANI, where the genome IDs match column 3 of "
			f"busco_taxid_genome_link.txt.\n"
			f"Pass --ani-file to point somewhere other than "
			f"{db_dir / 'ani.tsv'}."
		)
	elif Path(ani_path).stat().st_size == 0:
		raise ValueError(f"ANI file is empty: {ani_path}")


	ref_fasta = db_dir / f"{db_prefix}.fasta"
	if not ref_fasta.exists():
		raise ValueError(
			f"Reference marker FASTA not found: {ref_fasta}\n"
			f"RPKS needs this file's .fai index to compute each taxon's full "
			f"marker length; it cannot use the alignment or a precomputed "
			f"length file instead (see marker_lengths_from_reference in "
			f"eukfrac_calc.py for why)."
		)
	fai_path = Path(str(ref_fasta) + ".fai")
	if not fai_path.exists():
		raise ValueError(
			f"No .fai index for the reference FASTA: {fai_path}\n"
			f"Build it once with:\n"
			f"    samtools faidx {ref_fasta}\n"
			f"Rebuild it after any change to {ref_fasta} (added, removed, or "
			f"renamed sequences), or RPKS will silently use a stale length "
			f"for whatever changed."
		)
	elif fai_path.stat().st_size == 0:
		raise ValueError(f".fai index is empty: {fai_path}")

	logger.debug(f"Database validated: {db_dir}")


def _get_sample_paths(config: dict, sample: str, info) -> List[Path]:

	paired_end = config["paired_end"]
	fq_dir = Path(config["fq_dir"])

	if info and isinstance(info, dict) and "reads1" in info:
		# Absolute paths stored by ConfigBuilder — use them directly
		paths = [Path(info["reads1"])]
		if paired_end and "reads2" in info:
			paths.append(Path(info["reads2"]))
	else:
		# TSV / config-file mode: reconstruct from fq_dir + suffix
		if paired_end:
			paths = [
				fq_dir / f"{sample}{config['fwd_suffix']}",
				fq_dir / f"{sample}{config['rev_suffix']}",
			]
		else:
			paths = [fq_dir / f"{sample}{config['se_suffix']}"]

	return paths


def check_fastq_files(config: dict) -> None:
	samples = config["samples"]
	missing = []

	for sample, info in samples.items():
		for p in _get_sample_paths(config, sample, info):
			if not p.exists():
				missing.append(str(p))

	if missing:
		raise ValueError(
			"Missing input fastq files:\n  " + "\n  ".join(missing)
		)

	logger.debug(f"All fastq files found for {len(samples)} sample(s)")


def check_readlengths(config: dict) -> None:
	samples = config["samples"]
	expected_readlen = config["readlen"]
	warnings = []

	for sample, info in samples.items():
		for fastq_file in _get_sample_paths(config, sample, info):
			actual_readlen = _get_readlen(fastq_file)
			if actual_readlen == 0:
				continue
			if abs(actual_readlen - expected_readlen) > 10:
				warnings.append(
					f"{fastq_file.name}: actual={actual_readlen} bp, expected={expected_readlen} bp"
				)

	if warnings:
		logger.warning(
			"Read length mismatch detected:\n  " + "\n  ".join(warnings)
		)
		logger.warning(
			"Consider updating --readlen or checking your input files"
		)


def _get_readlen(fastq_path: Path, max_reads: int = 10000) -> int:	
	counter = 0
	bases = 0
	
	try:
		if str(fastq_path).endswith(".gz"):
			with gzip.open(fastq_path, "rt") as handle:
				for record in SeqIO.parse(handle, "fastq"):
					if counter >= max_reads:
						break
					counter += 1
					bases += len(record.seq)
		else:
			for record in SeqIO.parse(fastq_path, "fastq"):
				if counter >= max_reads:
					break
				counter += 1
				bases += len(record.seq)
		
		return int(bases / counter) if counter > 0 else 0

	except Exception as e:
		logger.warning(f"Could not read {fastq_path}: {e}")
		return 0


def check_alignment_outputs(config: dict, force: bool) -> None:
	output_dir = Path(config["output_dir"])
	samples = config["samples"].keys()
	
	existing = []
	for sample in samples:
		aln_file = output_dir / "aln" / f"{sample}_aln_q10_lenfilter.sorted.bam"
		if aln_file.exists():
			existing.append(str(aln_file))
	
	if existing and not force:
		raise ValueError(
			f"Alignment output files already exist. Use --force to overwrite:\n  " +
			"\n  ".join(existing)
		)
	
	if existing and force:
		logger.info(f"Removing {len(existing)} existing alignment file(s)")
		for f in existing:
			os.remove(f)


def check_filter_outputs(config: dict, force: bool) -> None:	
	output_dir = Path(config["output_dir"])
	samples = config["samples"].keys()
	
	existing = []
	for sample in samples:
		files = [
			output_dir / "filtering" / f"{sample}_aln_q10_lenfilter_complexityfilter_dupfilter.sorted.bam",
			output_dir / "filtering" / f"{sample}_aln_q10_lenfilter_complexityfilter_dupfilter.sorted.bam.bai",
			output_dir / "filtering" / f"{sample}_read_counts_and_mismatches.txt",
			output_dir / "filtering" / f"{sample}_all_hits_table.txt",
			output_dir / f"{sample}_filtered_hits_table.txt",
			output_dir / f"{sample}_filtered_hits_eukfrac.txt",
		]
		
		for f in files:
			if f.exists():
				existing.append(str(f))
	
	if existing and not force:
		raise ValueError(
			f"Filter output files already exist. Use --force to overwrite:\n  " +
			"\n  ".join(existing[:5]) + (f"\n  and {len(existing)-5} more" if len(existing) > 5 else "")
		)
	
	if existing and force:
		logger.info(f"Removing {len(existing)} existing filter file(s)")
		for f in existing:
			os.remove(f)


def check_alignment_inputs(config: dict) -> None:
	output_dir = Path(config["output_dir"])
	samples = config["samples"].keys()
	
	missing = []
	for sample in samples:
		aln_file = output_dir / "aln" / f"{sample}_aln_q10_lenfilter.sorted.bam"
		if not aln_file.exists():
			missing.append(str(aln_file))
	
	if missing:
		raise ValueError(
			f"Alignment files required for filter mode not found:\n  " +
			"\n  ".join(missing) +
			"\n\nRun alignment first with: eukdetect single --mode aln"
		)

