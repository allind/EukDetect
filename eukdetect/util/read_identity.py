#!/usr/bin/env python

import argparse
import logging
import sys

import pysam

logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO,
)
logger = logging.getLogger(__name__)


def alignment_identity(rec):
	cig = rec.cigartuples or []
	aligned = sum(l for op, l in cig if op in (0, 7, 8))      # M, =, X
	indel = sum(l for op, l in cig if op in (1, 2))           # I, D
	if aligned == 0:
		aligned = rec.query_alignment_length or 0
	if aligned == 0:
		return None, 0, 0
	try:
		nm = rec.get_tag("NM")
	except KeyError:
		return None, aligned, 0
	subs = max(nm - indel, 0)
	return (aligned - subs) / aligned, aligned, subs


def main(argv=None):
	p = argparse.ArgumentParser(description=__doc__,
								formatter_class=argparse.RawDescriptionHelpFormatter)
	p.add_argument("--bam", required=True,
				   help="The standard filtered, sorted BAM (q10, complexity, "
						"dedup) -- one alignment per read, as bowtie2 reports "
						"by default.")
	p.add_argument("--out", required=True, help="Per-read identity table")
	args = p.parse_args(argv)

	total = 0
	no_nm = 0
	unexpected_secondary = 0

	try:
		bam = pysam.AlignmentFile(args.bam, "rb")
	except Exception as e:
		logger.error(f"Could not open {args.bam}: {e}")
		sys.exit(1)

	with bam, open(args.out, "w") as out:
		out.write("Read\tMate\tMarker\tIdentity\tAligned_len\tMismatches\tAS\n")
		for rec in bam.fetch(until_eof=True):
			if rec.is_unmapped or rec.reference_name is None:
				continue
			if rec.is_secondary:
				unexpected_secondary += 1
				continue
			total += 1
			ident, aligned, subs = alignment_identity(rec)
			if ident is None:
				no_nm += 1
				continue
			try:
				score = rec.get_tag("AS")
			except KeyError:
				score = ""
			mate = 2 if rec.is_read2 else 1
			out.write(
				f"{rec.query_name}\t{mate}\t{rec.reference_name}\t"
				f"{ident * 100:.3f}\t{aligned}\t{subs}\t{score}\n"
			)

	if no_nm:
		logger.warning(
			f"{no_nm} of {total} alignments had no NM tag; identity could not "
			f"be computed for them and they are omitted."
		)
	if unexpected_secondary:
		logger.warning(
			f"{unexpected_secondary} alignment records were flagged secondary "
			f"in a BAM expected to hold one alignment per read; skipped."
		)
	logger.info(f"Wrote {total - no_nm} read alignments to {args.out}")


if __name__ == "__main__":
	main()
