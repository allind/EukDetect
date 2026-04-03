#!/usr/bin/env python

#complexity_filter.py - Filter low-complexity reads and remove duplicates from a BAM file.
#Replaces the bedtools bamtofastq | kz | remove_low_complexity | fixmate | markdup pipeline.
#Replaces kz from komplexity channel eclarke because not available in bioconda

#Complexity is measured as: unique 4-mers / sequence length  (identical to kz default).
#Reads with complexity < 0.5 are discarded

#After complexity filtering, PCR duplicates are removed using pysam's read-name-based


import sys
import os
import logging
import pysam
from collections import defaultdict

logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO,
)
logger = logging.getLogger(__name__)

COMPLEXITY_THRESHOLD = 0.5
KMER_SIZE = 4


def kmer_complexity(seq: bytes, k: int = KMER_SIZE) -> float:

	seq = seq.upper()
	n = len(seq)
	if n < k:
		return 0.0
	kmers = set()
	for i in range(n - k + 1):
		kmers.add(seq[i:i + k])
	return len(kmers) / n


def is_low_complexity(read: pysam.AlignedSegment, threshold: float = COMPLEXITY_THRESHOLD) -> bool:
	seq = read.query_sequence
	if seq is None:
		return False
	return kmer_complexity(seq.encode()) < threshold


def make_dup_key_paired(read: pysam.AlignedSegment):
	if read.is_read1:
		own = (read.reference_id, read.reference_start, read.is_reverse)
		mate = (read.next_reference_id, read.next_reference_start, read.mate_is_reverse)
	else:
		own = (read.next_reference_id, read.next_reference_start, read.mate_is_reverse)
		mate = (read.reference_id, read.reference_start, read.is_reverse)
	return (min(own, mate), max(own, mate))


def make_dup_key_single(read: pysam.AlignedSegment):
	return (read.reference_id, read.reference_start, read.is_reverse)


def filter_bam(input_bam: str, output_bam: str) -> None:

	with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inbam:
		header = inbam.header.to_dict()
		paired_end = False
		for read in inbam:
			if not read.is_unmapped:
				paired_end = read.is_paired
				break

	logger.info(f"Mode: {'paired-end' if paired_end else 'single-end'}")

	logger.info("Complexity filtering")

	if paired_end:
		pair_complexity: dict = defaultdict(list)
		with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inbam:
			for read in inbam:
				passes = not is_low_complexity(read)
				pair_complexity[read.query_name].append(passes)
		pass_names = {
			name for name, results in pair_complexity.items()
			if all(results)
		}
		logger.info(
			f"  {len(pass_names)} read pairs passed complexity filter "
			f"(discarded {len(pair_complexity) - len(pass_names)} pairs)"
		)
	else:
		pass_names = set()
		total = 0
		with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inbam:
			for read in inbam:
				total += 1
				if not is_low_complexity(read):
					pass_names.add(read.query_name)
		logger.info(
			f"  {len(pass_names)}/{total} reads passed complexity filter"
		)

	logger.info("Duplicate detection")

	best: dict = {}		   # key -> (query_name, quality_sum)
	dup_keys_seen: dict = defaultdict(list)

	with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inbam:
		for read in inbam:
			if read.query_name not in pass_names:
				continue
			if paired_end:
				key = make_dup_key_paired(read)
			else:
				key = make_dup_key_single(read)

			qual = read.query_qualities
			qsum = int(sum(qual)) if qual is not None else 0

			if key not in best or qsum > best[key][1]:
				best[key] = (read.query_name, qsum)
			dup_keys_seen[key].append(read.query_name)

	keep_names = {name for name, _ in best.values()}

	n_dup_positions = sum(
		1 for names in dup_keys_seen.values() if len(set(names)) > 1
	)
	logger.info(f"  Identified {n_dup_positions} duplicate positions")
	logger.info(f"  {len(keep_names)} reads/pairs kept after deduplication")

	logger.info("Writing filtered BAM")

	tmp_bam = output_bam + ".tmp_unsorted.bam"
	written = 0
	with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as inbam:
		with pysam.AlignmentFile(tmp_bam, "wb", header=header) as outbam:
			for read in inbam:
				if read.query_name in keep_names:
					outbam.write(read)
					written += 1

	logger.info(f"  Wrote {written} reads to temp BAM")

	logger.info("Sorting output BAM")
	pysam.sort("-o", output_bam, tmp_bam)
	os.unlink(tmp_bam)
	logger.info(f"Done. Output: {output_bam}")


def main():
	if len(sys.argv) != 3:
		print(f"Usage: {sys.argv[0]} <input.bam> <output.bam>", file=sys.stderr)
		sys.exit(1)

	input_bam = sys.argv[1]
	output_bam = sys.argv[2]

	logger.info(f"Input BAM:  {input_bam}")
	logger.info(f"Output BAM: {output_bam}")
	logger.info(f"Complexity threshold: < {COMPLEXITY_THRESHOLD} (kmer size k={KMER_SIZE})")

	filter_bam(input_bam, output_bam)


if __name__ == "__main__":
	main()
