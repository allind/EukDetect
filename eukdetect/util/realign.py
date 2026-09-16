#!/usr/bin/env python


import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

import pysam

try:
	from . import calibrate
except ImportError:
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	import calibrate

logger = logging.getLogger(__name__)

BUSCO_RE = re.compile(r"-\d+at\d+-")
DEFAULT_SCORE_TOLERANCE = 0


def marker_busco(seq_name):
	"""Extract the BUSCO id from a marker name, or a Collapse sentinel."""
	if "Collapse" in seq_name:
		return "Collapsed"
	m = BUSCO_RE.search(seq_name)
	return m.group(0).strip("-") if m else "Unknown"


def extract_cluster_reads(bam_path, marker_names, out_dir, prefix, paired_end):

	os.makedirs(out_dir, exist_ok=True)

	# qname -> {1: (seq, qual), 2: (seq, qual)}
	frags = defaultdict(dict)
	wanted = set(marker_names)

	with pysam.AlignmentFile(bam_path, "rb") as bam:
		present = wanted.intersection(bam.references)
		missing = wanted - present
		if missing:
			logger.debug(f"{len(missing)} cluster markers absent from BAM header; no reads there.")
		for ref in sorted(present):
			for read in bam.fetch(ref):
				if read.is_unmapped or read.query_sequence is None:
					continue
				mate = 2 if read.is_read2 else 1
				if mate in frags[read.query_name]:
					continue  # already captured this mate
				seq = read.query_sequence
				qual = read.qual if read.qual is not None else "I" * len(seq)

				if read.is_reverse:
					seq = _revcomp(seq)
					qual = qual[::-1]
				frags[read.query_name][mate] = (seq, qual)

	if not frags:
		return [], 0

	if paired_end:
		both = {q: m for q, m in frags.items() if 1 in m and 2 in m}
		orphan = {q: m for q, m in frags.items() if q not in both}
		r1 = os.path.join(out_dir, f"{prefix}_R1.fastq")
		r2 = os.path.join(out_dir, f"{prefix}_R2.fastq")
		with open(r1, "w") as f1, open(r2, "w") as f2:
			for q, m in both.items():
				f1.write(_fq(q, *m[1]))
				f2.write(_fq(q, *m[2]))
		paths = [r1, r2]
		if orphan:
			se = os.path.join(out_dir, f"{prefix}_SE.fastq")
			with open(se, "w") as f:
				for q, m in orphan.items():
					mate = m.get(1) or m.get(2)
					f.write(_fq(q, *mate))
			paths.append(se)
			logger.debug(f"{len(orphan)} half-mapped fragments written as single-end.")
		return paths, len(frags)

	se = os.path.join(out_dir, f"{prefix}_SE.fastq")
	with open(se, "w") as f:
		for q, m in frags.items():
			mate = m.get(1) or m.get(2)
			f.write(_fq(q, *mate))
	return [se], len(frags)


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s):
	return s.translate(_COMP)[::-1]


def _fq(name, seq, qual):
	return f"@{name}\n{seq}\n+\n{qual}\n"



def _index_key(marker_names):
	h = hashlib.sha1()
	for name in sorted(marker_names):
		h.update(name.encode())
		h.update(b"\n")
	return h.hexdigest()[:16]


def build_mini_index(marker_names, ref_fasta, cache_dir, cluster_id, bowtie2_build="bowtie2-build"):

	os.makedirs(cache_dir, exist_ok=True)
	key = _index_key(marker_names)
	prefix = os.path.join(cache_dir, f"{cluster_id}_{key}")
	done_marker = prefix + ".done"

	if os.path.exists(done_marker):
		logger.debug(f"Reusing cached mini index {prefix}")
		return prefix

	fasta_path = prefix + ".fasta"
	n_written = 0
	with pysam.FastaFile(ref_fasta) as ref, open(fasta_path, "w") as out:
		available = set(ref.references)
		for name in sorted(marker_names):
			if name not in available:
				logger.warning(f"Marker {name} not in {ref_fasta}; skipped.")
				continue
			out.write(f">{name}\n{ref.fetch(name)}\n")
			n_written += 1

	if n_written == 0:
		raise RuntimeError(f"No markers for {cluster_id} found in {ref_fasta}")

	cmd = [bowtie2_build, "--quiet", "--threads", "1", fasta_path, prefix]
	proc = subprocess.run(cmd, capture_output=True, text=True)
	if proc.returncode != 0:
		raise RuntimeError(
			f"bowtie2-build failed for {cluster_id}:\n{proc.stderr[:2000]}"
		)

	open(done_marker, "w").close()
	logger.debug(f"Built mini index {prefix} ({n_written} markers)")
	return prefix


def run_realignment(index_prefix, fastqs, out_sam, threads=1, bowtie2="bowtie2"):

	base = [
		bowtie2, "-a", "--end-to-end", "--very-sensitive",
		"--quiet", "--omit-sec-seq", "--no-unal",
		"-p", str(threads), "-x", index_prefix,
	]

	paired = [f for f in fastqs if f.endswith(("_R1.fastq", "_R2.fastq"))]
	singles = [f for f in fastqs if f.endswith("_SE.fastq")]

	sams = []
	tmp_idx = 0
	if len(paired) == 2:
		tmp = f"{out_sam}.part{tmp_idx}"
		cmd = base + ["--no-discordant", "-X", "1000", "-1", paired[0], "-2", paired[1], "-S", tmp]
		_run(cmd, "bowtie2 (paired realignment)")
		sams.append(tmp)
		tmp_idx += 1
	for se in singles:
		tmp = f"{out_sam}.part{tmp_idx}"
		cmd = base + ["-U", se, "-S", tmp]
		_run(cmd, "bowtie2 (single-end realignment)")
		sams.append(tmp)
		tmp_idx += 1

	if len(sams) == 1:
		shutil.move(sams[0], out_sam)
	else:
		with open(out_sam, "w") as out:
			for i, s in enumerate(sams):
				with open(s) as f:
					for line in f:
						if line.startswith("@") and i > 0:
							continue
						out.write(line)
				os.remove(s)
	return out_sam


def _run(cmd, what):
	proc = subprocess.run(cmd, capture_output=True, text=True)
	if proc.returncode != 0:
		raise RuntimeError(f"{what} failed:\n{' '.join(cmd)}\n{proc.stderr[:2000]}")


def score_fragments(sam_path, seq_taxid, score_tolerance=DEFAULT_SCORE_TOLERANCE):
	
	# qname -> marker -> mate -> (AS, matches, aligned_length)
	frag = defaultdict(lambda: defaultdict(dict))

	missing_nm = 0
	with pysam.AlignmentFile(sam_path, "r") as sam:
		for rec in sam:
			if rec.is_unmapped:
				continue
			try:
				score = rec.get_tag("AS")
			except KeyError:
				continue
			marker = rec.reference_name
			if marker not in seq_taxid:
				continue

			cig = rec.cigartuples or []
			alen = sum(l for op, l in cig if op in (0, 7, 8))   # M, =, X
			indel = sum(l for op, l in cig if op in (1, 2))     # I, D
			if alen == 0:
				alen = rec.query_alignment_length or 0
			try:
				nm = rec.get_tag("NM")
			except KeyError:
				nm = 0
				missing_nm += 1
			substitutions = max(nm - indel, 0)
			matches = max(alen - substitutions, 0)

			mate = 2 if rec.is_read2 else 1
			cur = frag[rec.query_name][marker].get(mate)
			if cur is None or matches > cur[1]:
				frag[rec.query_name][marker][mate] = (score, matches, alen)

	if missing_nm:
		logger.warning(
			f"{missing_nm} alignments lacked an NM tag; their identity was "
			f"treated as perfect, which makes the unique-evidence gate more "
			f"permissive rather than less."
		)

	results = {}
	for qname, markers in frag.items():
		species_score, species_ident, best_marker = {}, {}, {}
		for marker, mates in markers.items():
			taxid = seq_taxid[marker]
			matches = sum(v[1] for v in mates.values())
			alen = sum(v[2] for v in mates.values())
			if taxid not in species_score or matches > species_score[taxid]:
				species_score[taxid] = matches
				species_ident[taxid] = (matches / alen) if alen else 0.0
				best_marker[taxid] = marker
		if not species_score:
			continue
		top = max(species_score.values())
		winners = sorted(
			t for t, sc in species_score.items() if sc >= top - score_tolerance
		)
		results[qname] = {
			"winners": winners,
			"scores": species_score,
			"identity": species_ident,
			"best_marker": best_marker,
		}
	return results


def summarize_evidence(fragment_results, marker_busco_fn=None, **kw):

	ev, _ = calibrate.summarize_evidence(
		fragment_results, marker_busco_fn or marker_busco, **kw
	)
	return ev


def realign_cluster(
	cluster_id,
	marker_names,
	seq_taxid,
	bam_path,
	ref_fasta,
	cache_dir,
	paired_end,
	work_dir=None,
	threads=1,
	score_tolerance=DEFAULT_SCORE_TOLERANCE,
	anchor_identity=None,
	min_separation=None,
	anchor_quantile=None,
	anchor_slack=None,
	debug_dir=None,
):

	if debug_dir:
		work_dir = os.path.join(debug_dir, cluster_id)
		os.makedirs(work_dir, exist_ok=True)
		own_tmp = False
	else:
		own_tmp = work_dir is None
		work_dir = work_dir or tempfile.mkdtemp(prefix=f"eukdetect_{cluster_id}_")
	try:
		fastqs, n_frags = extract_cluster_reads(
			bam_path, marker_names, work_dir, cluster_id, paired_end
		)
		if n_frags == 0:
			logger.debug(f"{cluster_id}: no reads, skipping realignment.")
			return {}, {}, 0, {}

		index_prefix = build_mini_index(marker_names, ref_fasta, cache_dir, cluster_id)
		sam = os.path.join(work_dir, f"{cluster_id}.sam")
		run_realignment(index_prefix, fastqs, sam, threads=threads)

		if debug_dir:

			try:
				shutil.copy(index_prefix + ".fasta",
							os.path.join(work_dir, f"{cluster_id}_markers.fasta"))
			except Exception as e:
				logger.debug(f"Could not copy mini index FASTA: {e}")
			with open(os.path.join(work_dir, f"{cluster_id}_markers.txt"), "w") as f:
				for m in sorted(marker_names):
					f.write(f"{m}\t{seq_taxid.get(m, 'NA')}\n")

		frag_results = score_fragments(sam, seq_taxid, score_tolerance)
		kw = {}
		if anchor_identity is not None:
			kw["anchor_identity"] = anchor_identity
		if min_separation is not None:
			kw["min_separation"] = min_separation
		if anchor_quantile is not None:
			kw["anchor_quantile"] = anchor_quantile
		if anchor_slack is not None:
			kw["anchor_slack"] = anchor_slack
		evidence, diagnostics = calibrate.summarize_evidence(
			frag_results, marker_busco, cluster_id=cluster_id,
			debug_dir=work_dir if debug_dir else None, **kw
		)
		logger.info(
			f"{cluster_id}: realigned {n_frags} fragments across "
			f"{len(marker_names)} markers -> {len(evidence)} species with evidence"
		)
		diagnostics["n_fragments"] = n_frags
		diagnostics["n_markers"] = len(marker_names)
		return frag_results, evidence, n_frags, diagnostics
	finally:
		if own_tmp:
			shutil.rmtree(work_dir, ignore_errors=True)
