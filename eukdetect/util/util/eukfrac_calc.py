#!/usr/bin/env python
"""
eukfrac_calc.py - Taxonomic assignment and abundance estimation.

This replaces the previous two-pass heuristic. The old design chose a "primary"
taxon by read count and then transferred every read from any relative whose
pooled percent identity was no better. Two things were wrong with that. Pooled
PID is a property of a reference pileup rather than of a read, so it cannot say
where a read came from; and the transfer was all-or-nothing and ordered by
abundance, so a genuinely present minor species was absorbed by an abundant
relative and disappeared.

The pipeline is now:

  0. require    <database_dir>/ani.tsv, which decides which species compete
  1. load       database tables and per-marker read counts
  2. cluster    competition groups from ANI links UNIONED with NCBI
                rank-mates, so a wrong or over-split genus is corrected by
                ANI and a gap in ANI coverage is covered by taxonomy
                                                            (ani_clusters)
  3. nominate   any taxon with at least one read
  4. realign    competing candidates against a mini index   (realign)
  5. calibrate  measure spillover from bridge reads, then test whether each
                orphan read belongs where it landed or attaches to an anchor
                                                            (calibrate)
  6. assign     presence gate, EM, LCA fallback             (assign)
  6. report     the three output tables

Steps 4 and 5 only run for groups holding two or more candidates. A group with
one candidate has nothing to compete against and is reported directly, which is
the common case and costs nothing.

Output columns from the previous version are preserved so existing parsers keep
working, and new columns are appended. Two retained columns changed meaning:
Reads_aligned is now the reads only this taxon can explain, and Reads_reassigned
is now its EM share of reads that tied with a relative. Under the old code those
meant "reads that landed here" and "reads confiscated from a relative", which
are different quantities. Prefer the explicit Unique_reads, Shared_reads and
Resolution columns.
"""

import argparse
import logging
import os
import re
import sys
import textwrap
from collections import defaultdict

from ete3 import NCBITaxa

# This module is executed two ways: imported as part of the eukdetect package,
# and run directly as a script by the Snakemake rule (`python .../eukfrac_calc.py`).
# A bare relative import fails in the second case, so fall back to loading the
# sibling modules from this file's own directory.
try:
	from . import ani_clusters
	from . import realign as realign_mod
	from . import assign as assign_mod
except ImportError:
	sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
	import ani_clusters
	import realign as realign_mod
	import assign as assign_mod

logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO,
)
logger = logging.getLogger(__name__)

ORDERED_LABELS = ["phylum", "class", "order", "family", "genus", "species"]
BUSCO_RE = re.compile(r"-\d+at\d+-")


# ----------------------------------------------------------------------------
# 1. Loading
# ----------------------------------------------------------------------------

def load_gene_lengths(path):
	"""Cumulative marker length keyed by either taxid (numeric) or genome name."""
	taxid_genelen, genome_genelen = {}, {}
	try:
		with open(path) as f:
			for line in f:
				parts = line.strip().split("\t")
				if len(parts) < 2:
					continue
				ident = parts[0]
				try:
					length = int(parts[1])
				except ValueError:
					continue
				if ident.isdigit():
					taxid_genelen[ident] = length
				else:
					genome_genelen[ident] = length
	except FileNotFoundError:
		logger.error(f"Gene lengths file not found: {path}")
		sys.exit(1)
	logger.info(
		f"Loaded gene lengths for {len(taxid_genelen)} taxa and "
		f"{len(genome_genelen)} genomes"
	)
	return taxid_genelen, genome_genelen


def load_inherited_markers(path):
	"""taxid -> [buscos, specific_count, specific_and_inherited_count]"""
	full_seq_taxids = {}
	try:
		with open(path) as f:
			for line in f:
				parts = line.strip().split("\t")
				if len(parts) < 3:
					continue
				buscos = parts[1].split(",")
				full_seq_taxids[parts[0]] = [
					buscos, len(buscos), len(parts[2].split(","))
				]
	except FileNotFoundError:
		logger.error(f"Inherited markers file not found: {path}")
		sys.exit(1)
	logger.info(f"Loaded inherited markers for {len(full_seq_taxids)} taxa")
	return full_seq_taxids


def marker_busco(seq):
	if "Collapse" in seq:
		return "Collapsed"
	m = BUSCO_RE.search(seq)
	return m.group(0).strip("-") if m else "Unknown"


def parse_read_counts(path, seq_taxid):
	"""
	Read the per-marker count table.

	Returns (taxid_counts, n_records), where taxid_counts maps taxid to a list of
	[seq, count, correct_bases, total_bases, subjlen, coverage, pid, busco].

	Unlike the previous version this does no lineage lookup. Grouping is decided
	by ANI now, and resolving a genus for every marker line was a large part of
	the old runtime.
	"""
	taxid_counts = defaultdict(list)
	n = 0
	unknown_seqs = 0
	try:
		with open(path) as f:
			f.readline()  # header
			for line in f:
				n += 1
				parts = line.strip().split("\t")
				if len(parts) < 8:
					logger.warning(f"Skipping malformed line {n} in {path}")
					continue
				seq = parts[0]
				if seq not in seq_taxid:
					unknown_seqs += 1
					continue
				try:
					count = int(parts[1])
					correct_bases = int(parts[2])
					total_bases = int(parts[4])
					subjlen = int(parts[5])
					coverage = float(parts[6])
					pid = float(parts[7])
				except (ValueError, IndexError) as e:
					logger.warning(f"Error parsing line {n}: {e}; skipping")
					continue
				taxid_counts[seq_taxid[seq]].append(
					[seq, count, correct_bases, total_bases, subjlen,
					 coverage, pid, marker_busco(seq)]
				)
	except FileNotFoundError:
		logger.error(f"Read counts file not found: {path}")
		sys.exit(1)

	if unknown_seqs:
		logger.warning(
			f"{unknown_seqs} markers in the read counts file are absent from the "
			f"taxid link file and were skipped."
		)
	logger.info(f"Processed {n} alignment records over {len(taxid_counts)} taxa")
	return taxid_counts, n


def compute_taxon_stats(taxid_counts, taxid_seqs, seq_genomes,
						taxid_genelen, genome_genelen, min_identity_pct=97.0):
	"""
	Per-taxon aggregates from the raw alignment, before any competition.

	These fill the all-hits table and serve as the fallback for taxa that never
	get realigned.
	"""
	stats, observed_genomes = {}, defaultdict(set)
	for tax, entries in taxid_counts.items():
		reads = sum(e[1] for e in entries)
		correct = sum(e[2] for e in entries)
		total_bases = sum(e[3] for e in entries)
		subj_len = sum(e[4] for e in entries)
		buscos = {e[7] for e in entries} - {"Collapsed", "Unknown"}

		for e in entries:
			for g in seq_genomes.get(e[0], ()):
				observed_genomes[tax].add(g)

		# Prefer a genome-specific marker length when exactly one genome was
		# seen, since that is the length the reads could actually have come
		# from. Otherwise fall back to the species-level total.
		genomes = observed_genomes[tax]
		marker_length = 0
		if len(genomes) == 1:
			marker_length = genome_genelen.get(next(iter(genomes)), 0)
		if not marker_length:
			marker_length = taxid_genelen.get(tax, subj_len)

		# Evidence restricted to markers whose reads actually match this taxon
		# well. A taxon that never has to compete still needs an absolute
		# quality bar, otherwise a species with no close relative in the
		# database is called on a handful of 95%-identity strays from something
		# else entirely. Counted per marker rather than on the pooled identity,
		# so a few junk markers cannot drag down a genuinely present taxon.
		hi_buscos, hi_reads = set(), 0
		hi_correct, hi_total = 0, 0
		for e in entries:
			if e[6] >= min_identity_pct and e[7] not in ("Collapsed", "Unknown"):
				hi_buscos.add(e[7])
				hi_reads += e[1]
				hi_correct += e[2]
				hi_total += e[3]

		stats[tax] = {
			"observed_markers": len(entries),
			"reads": reads,
			"buscos": buscos,
			"hi_identity_buscos": hi_buscos,
			"hi_identity_reads": hi_reads,
			# Identity of the qualifying markers alone. Reporting the pooled
			# figure for a species-level call would understate it whenever a few
			# repetitive markers drag the average down.
			"hi_identity_pid": round(hi_correct / hi_total * 100, 2) if hi_total else 0,
			"percent_identity": round(correct / total_bases * 100, 2) if total_bases else 0,
			"coverage": round(total_bases / subj_len * 100, 2) if subj_len else 0,
			"marker_length": marker_length,
			"total_markers": len(taxid_seqs.get(tax, [])),
		}
	return stats, observed_genomes


# ----------------------------------------------------------------------------
# 2-5. Competition
# ----------------------------------------------------------------------------

def resolve_all_clusters(clusters, taxid_cluster, provenance, taxid_counts,
						 taxon_stats, taxid_seqs, taxid_genelen, ncbi, args):
	"""
	Run competition for every group holding more than one candidate.

	The mini index is built from the markers of *candidate* species only, not
	from every member of the group. Group size reflects how many relatives the
	database happens to contain, whereas realignment cost should reflect what is
	actually in the sample: a 300-member genus with three species present should
	cost the same as a three-member group.

	Returns (calls_by_taxid, realigned_cluster_ids).
	"""
	candidates_by_cluster = defaultdict(list)
	for tax in taxid_counts:
		candidates_by_cluster[taxid_cluster.get(tax)].append(tax)

	multi = {c: t for c, t in candidates_by_cluster.items()
			 if c is not None and len(t) > 1}
	logger.info(
		f"{len(candidates_by_cluster)} groups have candidates; "
		f"{len(multi)} hold competing candidates and will be realigned"
	)

	calls, realigned, reports = {}, set(), []
	for cid, cand in sorted(multi.items()):
		marker_names, seq_taxid_local = [], {}
		for tax in cand:
			for s in taxid_seqs.get(tax, []):
				marker_names.append(s)
				seq_taxid_local[s] = tax
		marker_names = sorted(set(marker_names))
		if not marker_names:
			continue

		try:
			frag_results, evidence, n_frags, diag = realign_mod.realign_cluster(
				cluster_id=cid,
				marker_names=marker_names,
				seq_taxid=seq_taxid_local,
				bam_path=args.prefilter_bam,
				ref_fasta=args.reference_fasta,
				cache_dir=args.index_cache,
				paired_end=args.paired_end,
				threads=args.threads,
				score_tolerance=args.score_tolerance,
				anchor_identity=args.anchor_identity,
				min_separation=args.min_separation,
			)
		except Exception as e:
			# One failing group must not cost the whole sample. Leave its
			# candidates to the uncompeted path and say so loudly.
			logger.error(
				f"Realignment failed for {cid} ({','.join(cand)}): {e}. "
				f"Its candidates fall back to uncompeted reporting."
			)
			continue

		if not evidence:
			continue
		diag["cluster_id"] = cid
		diag["candidates"] = sorted(cand)
		diag["evidence"] = {
			t: {k: (sorted(v) if isinstance(v, set) else v)
				for k, v in e.items() if k != "identity_sum"}
			for t, e in evidence.items()
		}
		reports.append(diag)

		marker_lengths = {
			t: taxid_genelen.get(t, taxon_stats[t]["marker_length"]) for t in cand
		}
		for call in assign_mod.resolve_cluster(
			cluster_id=cid,
			member_taxids=cand,
			fragment_results=frag_results,
			evidence=evidence,
			marker_lengths=marker_lengths,
			ncbi=ncbi,
			min_unique_reads=args.min_unique_reads,
			min_unique_markers=args.min_unique_markers,
			# EM must count exactly the fragments that became evidence. Anything
			# below the anchor identity was already credited as attached, so
			# letting EM count it as well double-reports it.
			min_identity=args.anchor_identity,
		):
			call["provenance"] = provenance.get(cid, "ani")
			calls[call["taxid"]] = call
		realigned.add(cid)

	return calls, realigned, reports


def uncompeted_call(tax, cid, prov, stats, args):
	"""
	Build a call for a taxon that never had to compete.

	Two outcomes, matching what the competed path does:

	  no_competition     enough evidence at within-species identity, so the
	                     taxon is identified to species
	  nearest_reference  enough evidence overall but not at within-species
	                     identity, so something is present and this is the
	                     closest thing the database has. Reported with the
	                     identity its reads actually achieved rather than
	                     dropped, since dropping it loses a real organism.

	The reported identity is taken from whichever evidence supports the call:
	the qualifying markers for a species-level call, the pooled figure for a
	nearest-reference one. Both are read-weighted, so the number in the table is
	what the reads really matched at.
	"""
	base = {
		"taxid": tax,
		"cluster_id": cid or "NA",
		"provenance": prov,
		"shared_reads": 0.0,
		"attached_reads": 0,
		"spillover_reads": 0,
		"undecidable_reads": 0,
		"delta_to_anchor": None,
		"bridge_reads": 0,
		"members": [tax],
	}

	if (len(stats["hi_identity_buscos"]) >= args.min_unique_markers
			and stats["hi_identity_reads"] >= args.min_unique_reads):
		base.update({
			"resolution": "no_competition",
			"is_anchor": True,
			"unique_reads": stats["hi_identity_reads"],
			"unique_markers": len(stats["hi_identity_buscos"]),
			"mean_identity": stats["hi_identity_pid"],
			"assigned_reads": float(stats["hi_identity_reads"]),
		})
		return base

	if (len(stats["buscos"]) >= args.min_unique_markers
			and stats["reads"] >= args.min_unique_reads):
		base.update({
			"resolution": "nearest_reference",
			"is_anchor": False,
			"unique_reads": stats["reads"],
			"unique_markers": len(stats["buscos"]),
			"mean_identity": stats["percent_identity"],
			"assigned_reads": float(stats["reads"]),
		})
		return base

	return None


# ----------------------------------------------------------------------------
# 6. Reporting
# ----------------------------------------------------------------------------

def build_lineages(ncbi, taxids):
	"""Rank-labelled lineage strings, in the format the previous version emitted."""
	out = {}
	for tax in taxids:
		try:
			lineage_list = ncbi.get_lineage(int(tax))
			names = ncbi.get_taxid_translator(lineage_list)
			ranks = ncbi.get_rank(lineage_list)
			ranks_rev = {v: k for k, v in ranks.items()}
			lin = ""
			for label in ORDERED_LABELS:
				if label in ranks_rev:
					lin += f"{label}-{names[ranks_rev[label]]}|"
			out[str(tax)] = lin.strip("|").replace(" ", "_")
		except Exception as e:
			logger.warning(f"Could not build lineage for {tax}: {e}")
	return out


def rank_and_name(ncbi, tax):
	try:
		rank = list(ncbi.get_rank([tax]).values())[0]
		name = list(ncbi.get_taxid_translator([tax]).values())[0]
	except Exception:
		return None, None
	if rank == "no rank":
		try:
			parent = ncbi.get_lineage(tax)[-2]
			rank = list(ncbi.get_rank([parent]).values())[0]
		except Exception:
			rank = "no rank"
	return rank, name


def write_alltab(path, ncbi, taxon_stats, observed_genomes, passing,
				 taxid_cluster, provenance):
	order = sorted(taxon_stats, key=lambda t: taxon_stats[t]["observed_markers"],
				   reverse=True)
	with open(path, "w") as dest:
		dest.write(
			"Name\tTaxid\tRank\tObserved_markers\tRead_counts\t"
			"Total_marker_coverage\tPercent_identity\tTotal_marker_length\t"
			"Genomes\tFiltered\tCluster_id\tCluster_provenance\n"
		)
		for tax in order:
			rank, name = rank_and_name(ncbi, tax)
			if rank is None:
				continue
			s = taxon_stats[tax]
			genomes = sorted(observed_genomes.get(tax, set()))
			cid = taxid_cluster.get(tax, "NA")
			dest.write(
				f"{name}\t{tax}\t{rank}\t{s['observed_markers']}\t{s['reads']}\t"
				f"{s['coverage']}%\t{s['percent_identity']}%\t{s['marker_length']}\t"
				f"{','.join(genomes) if genomes else 'NA'}\t"
				f"{'No' if tax in passing else 'Yes'}\t"
				f"{cid}\t{provenance.get(cid, 'NA')}\n"
			)
	logger.info(f"Wrote all-hits table to {path}")


def call_marker_length(tax, call, taxid_genelen, taxon_stats):
	"""
	Marker length for one call.

	An unresolved-cluster call sits at an internal rank that has no entry in the
	gene-length table, so it carries its own length hint from the members it
	stands in for.
	"""
	if call.get("marker_length"):
		return call["marker_length"]
	return taxid_genelen.get(
		tax, taxon_stats.get(tax, {}).get("marker_length", 0)
	)


def compute_abundances(ncbi, calls, taxid_genelen, taxon_stats, observed_genomes):
	"""
	RPKS and relative abundance over the tree of called taxa.

	The unit of abundance is a *call*, not a taxonomic rank. The previous version
	normalized over "the lowest rank that has data", which works while every call
	is a species but silently zeroes out any call made at an internal rank. An
	unresolved-cluster call is exactly that: it represents one real organism and
	its reads must count, even though it is reported at a genus or family node.

	Internal nodes sum their own call, if any, plus everything below them, so a
	genus holding both a resolved species and an unresolved sibling group totals
	correctly rather than double counting.
	"""
	passing = sorted(calls)
	try:
		# Build the node set from the full lineages rather than calling
		# get_topology on the taxids alone. ete3 takes a separate code path when
		# given exactly one taxid that requires a <db>.traverse.pkl cache, which
		# only exists if the database was built by ete3's own updater. A sample
		# with a single detected species is entirely normal, so relying on that
		# path makes a common case crash. Passing the lineages guarantees more
		# than one node and also saves building the tree twice.
		all_ids = set()
		for t in passing:
			try:
				all_ids.update(ncbi.get_lineage(int(t)))
			except Exception as e:
				logger.warning(f"No lineage for called taxon {t}: {e}")
		if not all_ids:
			logger.error("No lineages could be resolved for any called taxon.")
			sys.exit(1)
		tree = ncbi.get_topology(sorted(all_ids), intermediate_nodes=True)
	except SystemExit:
		raise
	except Exception as e:
		logger.error(f"Error building tree: {e}")
		sys.exit(1)

	# Per-call RPKS, then relative abundance across all calls.
	unit_rpks, unit_len = {}, {}
	for tax, call in calls.items():
		L = call_marker_length(tax, call, taxid_genelen, taxon_stats)
		unit_len[tax] = L
		unit_rpks[tax] = (call["assigned_reads"] / (L / 1000)) if L > 0 else 0.0
		if L <= 0:
			logger.warning(
				f"Taxon {tax} has no marker length; its "
				f"{call['assigned_reads']} reads cannot contribute to relative "
				f"abundance."
			)
	total_rpks = sum(unit_rpks.values())

	node_reads = defaultdict(float)
	node_len = defaultdict(float)
	node_rel = defaultdict(float)
	node_genomes = defaultdict(set)

	for node in tree.traverse("postorder"):
		name = node.name
		if name in calls:
			node_reads[name] += calls[name]["assigned_reads"]
			node_len[name] += unit_len[name]
			node_rel[name] += (
				(unit_rpks[name] / total_rpks * 100) if total_rpks > 0 else 0.0
			)
			node_genomes[name] |= observed_genomes.get(name, set())
			# An unresolved-cluster call sits on a taxid that never had reads of
			# its own, so its genomes come from the members it stands in for.
			for member in calls[name].get("members", ()):
				node_genomes[name] |= observed_genomes.get(member, set())
		for child in node.children:
			node_reads[name] += node_reads[child.name]
			node_len[name] += node_len[child.name]
			node_rel[name] += node_rel[child.name]
			node_genomes[name] |= node_genomes[child.name]

	relabs = {}
	for node in tree.traverse():
		name = node.name
		if node_reads.get(name, 0) <= 0 and name not in calls:
			continue
		L = node_len.get(name, 0)
		relabs[name] = [
			(node_reads[name] / (L / 1000)) if L > 0 else 0.0,
			L,
			node_rel.get(name, 0.0),
			node_reads.get(name, 0.0),
		]

	return tree, relabs, node_genomes


def write_primarytab(path, ncbi, calls, lineages, taxon_stats, observed_genomes,
					 taxid_genelen):
	with open(path, "w") as dest:
		dest.write(
			"Name\tRank\tLineage\tTaxid\tTotal_reads\tTotal_marker_length\t"
			"RPKS\tReads_aligned\tPID_aligned\tGenomes\t"
			"Reads_reassigned\tPID_reassigned\tReassigned_genomes\t"
			"Cluster_id\tCluster_provenance\tResolution\t"
			"Unique_reads\tUnique_markers\tShared_reads\t"
			"Is_anchor\tRead_identity\tDelta_to_anchor\tBridge_reads\t"
			"Attached_reads\tSpillover_reads\tUndecidable_reads\t"
			"Competing_taxa\n"
		)
		for tax in sorted(calls, key=lambda t: -calls[t]["assigned_reads"]):
			call = calls[tax]
			rank, name = rank_and_name(ncbi, tax)
			if rank is None:
				continue
			s = taxon_stats.get(tax, {})
			marker_len = call_marker_length(tax, call, taxid_genelen, taxon_stats)
			total_reads = call["assigned_reads"]
			rpks = round(total_reads / (marker_len / 1000), 4) if marker_len else 0.0
			genomes = sorted(observed_genomes.get(tax, set()))
			competitors = [m for m in call["members"] if m != tax]
			# Everything beyond what is uniquely this taxon's, i.e. its EM share
			# of the tied reads. Kept under the old column name.
			em_share = round(max(total_reads - call["unique_reads"], 0.0), 3)
			# PID is a pileup statistic on the original alignment and has no
			# meaningful per-competitor version, so none is invented here.
			dest.write(
				f"{name}\t{rank}\t{lineages.get(tax, rank)}\t{tax}\t"
				f"{round(total_reads, 2)}\t{int(marker_len)}\t{rpks}\t"
				f"{call['unique_reads']}\t{s.get('percent_identity', 0)}%\t"
				f"{','.join(genomes) if genomes else 'NA'}\t"
				f"{em_share}\tNA\t"
				f"{','.join(competitors) if competitors else 'None'}\t"
				f"{call['cluster_id']}\t{call.get('provenance', 'NA')}\t"
				f"{call['resolution']}\t"
				f"{call['unique_reads']}\t{call['unique_markers']}\t"
				f"{call['shared_reads']}\t"
				f"{'yes' if call.get('is_anchor') else 'no'}\t"
				f"{call.get('mean_identity', 0.0)}%\t"
				f"{call.get('delta_to_anchor') if call.get('delta_to_anchor') is not None else 'NA'}\t"
				f"{call.get('bridge_reads', 0)}\t"
				f"{call.get('attached_reads', 0)}\t"
				f"{call.get('spillover_reads', 0)}\t"
				f"{call.get('undecidable_reads', 0)}\t"
				f"{','.join(competitors) if competitors else 'None'}\n"
			)
	logger.info(f"Wrote primary table to {path}")


def write_eukfrac(path, ncbi, tree, relabs, lineages, calls, node_genomes):
	with open(path, "w") as dest:
		dest.write(
			"Lineage\tRank\tName\tTaxID\tRPKS\tRelative_abundance\t"
			"Reads_total\tTotal_marker_length\tGenomes\tResolution\n"
		)
		for node in tree.traverse("preorder"):
			rank, name = rank_and_name(ncbi, node.name)
			if rank is None:
				continue
			if rank == "no rank" and node.is_leaf():
				continue
			if node.name not in relabs:
				continue
			# Report RPKS for any node that is itself a call. Restricting this to
			# species rank would blank out unresolved-cluster calls, which are
			# real abundance units reported at an internal rank.
			rpks = round(relabs[node.name][0], 4) if node.name in calls else "NA"
			genomes = sorted(node_genomes.get(node.name, set()))
			res = calls.get(node.name, {}).get("resolution", "aggregate")
			lin = lineages.get(node.name)
			if lin is None:
				lin = build_lineages(ncbi, [node.name]).get(node.name, rank)
			dest.write(
				f"{lin}\t{rank}\t{name}\t{node.name}\t"
				f"{rpks}\t{round(relabs[node.name][2], 4)}\t"
				f"{round(relabs[node.name][3], 2)}\t{int(relabs[node.name][1])}\t"
				f"{','.join(genomes) if genomes else 'NA'}\t{res}\n"
			)
	logger.info(f"Wrote relative abundance table to {path}")


def write_realignment_report(path, reports, ncbi):
	"""
	Per-cluster record of what the realignment and calibration actually did.

	Without this the whole step is invisible: the mini index, the SAM and every
	per-fragment decision live in a temporary directory that is deleted when the
	cluster finishes. The two record types are `species`, one row per candidate
	with its identity distribution and where its fragments ended up, and `pair`,
	one row per anchor/candidate combination with the measured divergence and
	the number of bridge fragments it was measured from. A `pair` row with zero
	bridge fragments is the case where divergence could not be measured at all.
	"""
	with open(path, "w") as dest:
		dest.write(
			"Record\tCluster_id\tTaxid\tName\tIs_anchor\tAnchor_level\t"
			"Fragments\tAt_anchor_identity\tIdentity_p10\tIdentity_p50\t"
			"Identity_p90\tUnique\tAttached\tSpillover\tUndecidable\t"
			"Shared\tMarkers\tDelta\tBridge_fragments\tNote\n"
		)
		for rep in reports:
			cid = rep.get("cluster_id", "NA")
			anchors = rep.get("anchors", {})
			merged = rep.get("merged_anchors", {})
			prof = rep.get("identity_profile", {})
			ev = rep.get("evidence", {})

			for taxid in sorted(set(list(prof) + list(ev))):
				_, name = rank_and_name(ncbi, taxid)
				p = prof.get(taxid, {})
				e = ev.get(taxid, {})
				note = ""
				if taxid in merged:
					note = f"indistinguishable from {merged[taxid]}"
				dest.write(
					f"species\t{cid}\t{taxid}\t{name or 'NA'}\t"
					f"{'yes' if taxid in anchors else 'no'}\t"
					f"{anchors.get(taxid, 'NA')}\t"
					f"{p.get('n', 0)}\t{p.get('n_at_anchor_identity', 0)}\t"
					f"{p.get('p10', 'NA')}\t{p.get('p50', 'NA')}\t{p.get('p90', 'NA')}\t"
					f"{e.get('unique_reads', 0)}\t{e.get('attached_reads', 0)}\t"
					f"{e.get('spillover_reads', 0)}\t{e.get('undecidable_reads', 0)}\t"
					f"{round(e.get('shared_reads', 0.0), 3)}\t"
					f"{len(e.get('unique_buscos', []))}\t\t\t{note}\n"
				)

			for key, (delta, nbridge) in sorted(rep.get("deltas", {}).items()):
				anchor, cand = key.split("->", 1)
				note = "no bridge fragments" if not nbridge else ""
				dest.write(
					f"pair\t{cid}\t{cand}\tvs {anchor}\t\t\t\t\t\t\t\t"
					f"\t\t\t\t\t\t{delta if delta is not None else 'NA'}\t"
					f"{nbridge}\t{note}\n"
				)
	logger.info(f"Wrote realignment report to {path}")


def write_message(paths, message):
	for p in paths:
		try:
			with open(p, "w") as f:
				f.write(message + "\n")
		except Exception as e:
			logger.error(f"Error writing to {p}: {e}")


# ----------------------------------------------------------------------------
# Orchestration
# ----------------------------------------------------------------------------

def build_parser():
	p = argparse.ArgumentParser(
		description=textwrap.dedent("""\
			Assign reads to taxa and estimate abundance.

			Competing species are determined by genome ANI rather than NCBI
			genus, and resolved by realigning their reads against a mini
			database so presence is decided on read-level evidence.
		"""),
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	p.add_argument("--dbfile", required=True, help="ete3 NCBI taxonomy sqlite")
	p.add_argument("--inherited_markers", required=True)
	p.add_argument("--taxid_link", required=True)
	p.add_argument("--readcounts", required=True)
	p.add_argument("--taxid_genelens", required=True)
	p.add_argument("--eukfrac", required=True)
	p.add_argument("--primarytab", required=True)
	p.add_argument("--alltab", required=True)

	p.add_argument("--prefilter_bam", required=True,
				   help="Pre-MAPQ BAM used as the read source for realignment")
	p.add_argument("--reference_fasta", required=True,
				   help="Marker FASTA used to build mini indexes")
	p.add_argument("--ani_file", required=True,
				   help="genome_1<TAB>genome_2<TAB>ANI. A required database "
						"file, expected at <database_dir>/ani.tsv")
	p.add_argument("--ani_cutoff", type=float, default=0.80)
	p.add_argument("--ani_linkage", default="single",
				   choices=("single", "complete", "average"))
	p.add_argument("--group_by_taxonomy", action="store_true", default=True,
				   help="Union NCBI rank-mates into competition groups "
						"alongside ANI links (default on)")
	p.add_argument("--no_taxonomy_grouping", dest="group_by_taxonomy",
				   action="store_false")
	p.add_argument("--taxonomy_rank", default="genus",
				   help="NCBI rank unioned with ANI (default genus)")
	p.add_argument("--realign", action="store_true", default=True)
	p.add_argument("--no-realign", dest="realign", action="store_false")
	p.add_argument("--score_tolerance", type=int, default=0)
	p.add_argument("--index_cache", default="")
	p.add_argument("--min_unique_reads", type=int, default=4)
	p.add_argument("--min_unique_markers", type=int, default=2)
	p.add_argument("--anchor_identity", type=float, default=0.99,
				   help="Identity at which a species' own reads are taken to "
						"match it exactly, making it an anchor. A property of "
						"read length and error rate, not of taxonomy.")
	p.add_argument("--realign_report", default=None,
				   help="TSV recording what the realignment and calibration did "
						"for each cluster; without it the step leaves no trace.")
	p.add_argument("--min_separation", type=float, default=0.01,
				   help="Two species whose measured divergence is below this "
						"cannot be told apart by their reads; such fragments "
						"are reported as undecidable rather than guessed.")
	p.add_argument("--paired_end", default="true")
	p.add_argument("--threads", type=int, default=1)
	return p


def main(argv=None):
	args = build_parser().parse_args(None if argv is None else argv[1:])
	args.paired_end = str(args.paired_end).lower() in ("true", "1", "yes")
	if not args.index_cache:
		args.index_cache = os.path.join(
			os.path.dirname(os.path.abspath(args.alltab)), "realign_index_cache"
		)

	try:
		ncbi = NCBITaxa(args.dbfile)
	except Exception as e:
		logger.error(f"Failed to load NCBI taxonomy from {args.dbfile}: {e}")
		sys.exit(1)

	# --- 1. load -----------------------------------------------------------
	taxid_genelen, genome_genelen = load_gene_lengths(args.taxid_genelens)
	link = ani_clusters.load_genome_taxid_map(args.taxid_link)
	seq_taxid, seq_genomes = link["seq_taxid"], link["seq_genomes"]
	taxid_seqs = defaultdict(list)
	for s, t in seq_taxid.items():
		taxid_seqs[t].append(s)
	load_inherited_markers(args.inherited_markers)

	taxid_counts, n_records = parse_read_counts(args.readcounts, seq_taxid)
	if n_records == 0 or not taxid_counts:
		msg = "Empty read count file. Likely no aligned reads in sample."
		logger.warning(msg)
		write_message([args.eukfrac, args.alltab, args.primarytab], msg)
		if args.realign_report:
			write_realignment_report(args.realign_report, [], ncbi)
		return

	taxon_stats, observed_genomes = compute_taxon_stats(
		taxid_counts, taxid_seqs, seq_genomes, taxid_genelen, genome_genelen,
		min_identity_pct=args.anchor_identity * 100.0,
	)

	# --- 2. cluster --------------------------------------------------------
	# The ANI table is a required database file. Running without it would
	# complete and emit plausible-looking tables in which no species were ever
	# competed against each other, which is a far worse outcome than stopping.
	if not os.path.exists(args.ani_file):
		logger.error(
			f"ANI file not found: {args.ani_file}\n"
			f"It is a required database file, expected at "
			f"<database_dir>/{ani_clusters.ANI_FILENAME}, with three "
			f"tab-separated columns: genome_1, genome_2, ANI."
		)
		sys.exit(1)

	res = ani_clusters.build(
		args.ani_file, args.taxid_link,
		cutoff=args.ani_cutoff, linkage=args.ani_linkage,
		ncbi=ncbi if args.group_by_taxonomy else None,
		taxonomy_rank=args.taxonomy_rank,
		group_by_taxonomy=args.group_by_taxonomy,
		# Rank lineages come from <database_dir>/identity_cache.json when it is
		# present. The database directory is wherever the taxid link file lives.
		db_dir=os.path.dirname(os.path.abspath(args.taxid_link)),
	)
	clusters = res["clusters"]
	taxid_cluster = res["taxid_cluster"]
	provenance = res["provenance"]

	if not clusters:
		logger.error(
			f"No competition groups were built from {args.ani_file}. The genome "
			f"IDs in it probably do not match column 3 of {args.taxid_link}."
		)
		sys.exit(1)

	# --- 3-5. nominate, realign, assign ------------------------------------
	calls, realigned_clusters, realign_reports = {}, set(), []
	can_realign = bool(
		args.realign
		and os.path.exists(args.prefilter_bam)
		and os.path.exists(args.reference_fasta)
	)
	if args.realign and not can_realign:
		missing = [p for p in (args.prefilter_bam, args.reference_fasta)
				   if not os.path.exists(p)]
		logger.error(
			f"Realignment is enabled but required inputs are missing: "
			f"{', '.join(missing)}. Species cannot be disambiguated without "
			f"them; pass --no-realign to run without disambiguation."
		)
		sys.exit(1)

	if can_realign:
		calls, realigned_clusters, realign_reports = resolve_all_clusters(
			clusters, taxid_cluster, provenance, taxid_counts, taxon_stats,
			taxid_seqs, taxid_genelen, ncbi, args
		)

	# Anything competition did not settle is reported on its own evidence.
	# Candidates in a group that WAS realigned are skipped: they competed and
	# had no unique evidence, so re-admitting them here would undo the result.
	for tax in taxid_counts:
		if tax in calls or taxid_cluster.get(tax) in realigned_clusters:
			continue
		cid = taxid_cluster.get(tax)
		call = uncompeted_call(
			tax, cid, provenance.get(cid, "NA"), taxon_stats[tax], args
		)
		if call:
			calls[tax] = call

	if not calls:
		msg = "No taxa passing filter requirements."
		logger.warning(msg)
		write_alltab(args.alltab, ncbi, taxon_stats, observed_genomes,
					 set(), taxid_cluster, provenance)
		write_message([args.primarytab, args.eukfrac], msg)
		if args.realign_report:
			write_realignment_report(args.realign_report, realign_reports, ncbi)
		return

	n_species = sum(1 for c in calls.values() if c["resolution"] == "species")
	n_unres = sum(1 for c in calls.values() if c["resolution"] == "unresolved_cluster")
	logger.info(
		f"{len(calls)} taxa called: {n_species} resolved by competition, "
		f"{n_unres} unresolved groups reported at their LCA, "
		f"{len(calls) - n_species - n_unres} uncompeted"
	)

	# --- 6. report ---------------------------------------------------------
	write_alltab(args.alltab, ncbi, taxon_stats, observed_genomes,
				 set(calls), taxid_cluster, provenance)
	if args.realign_report:
		write_realignment_report(args.realign_report, realign_reports, ncbi)

	lineages = build_lineages(ncbi, list(calls))
	tree, relabs, node_genomes = compute_abundances(
		ncbi, calls, taxid_genelen, taxon_stats, observed_genomes
	)
	write_primarytab(args.primarytab, ncbi, calls, lineages, taxon_stats,
					 observed_genomes, taxid_genelen)
	write_eukfrac(args.eukfrac, ncbi, tree, relabs, lineages, calls, node_genomes)
	logger.info("Analysis complete.")


if __name__ == "__main__":
	main(sys.argv)
