#!/usr/bin/env python
"""
ani_clusters.py - Build competition groups from genome-genome ANI.

Replaces NCBI genus as the unit of taxonomic disambiguation. NCBI genus is a
nomenclatural label and does not track genome similarity: closely related
species are routinely split across genera, and divergent species are routinely
lumped into one. Using genus to decide which taxa compete for reads therefore
both misses real competitors and forces unrelated taxa to compete.

An "ANI cluster" is a set of species whose representative genomes are linked at
or above an ANI cutoff. Only species inside the same cluster ever compete for
reads. Species in no cluster are singletons and are reported as-is.

Input ANI file is 3 columns, no header:
    genome_1 <TAB> genome_2 <TAB> ANI

ANI may be on a 0-1 or 0-100 scale; the scale is detected and normalized to
0-1. The cutoff is always expressed on the 0-1 scale (default 0.80).

Linkage:
    single   (default) transitive closure over edges >= cutoff. Inclusive:
             it is cheap to let an extra species into a competition group,
             because the read-level realignment step resolves it. It is
             expensive to leave a true competitor out, because that is exactly
             the case where a low-abundance species gets absorbed by an
             abundant relative and disappears.
    complete every pair in the group must be >= cutoff. Conservative; use when
             single linkage is chaining unrelated genomes together.
    average  mean similarity between merged groups must be >= cutoff.

Species-level closure: a species may have several representative genomes. If
those genomes land in different genome clusters, the clusters are merged, since
a species must compete as a single unit.
"""

import argparse
import hashlib
import json
import logging
import os
import sys
from collections import defaultdict

logger = logging.getLogger(__name__)

# The ANI table is a required database file, expected at this name inside the
# database directory. Defined here so the pipeline, the validator and the
# taxonomy step all agree on one spelling.
ANI_FILENAME = "ani.tsv"

# Precomputed per-taxon rank lineages, built once by
# build_db/make_identity_cache.py and stored beside the other database files.
IDENTITY_CACHE_FILENAME = "identity_cache.json"

# 0.80, not the 0.95-0.97 usually quoted for species delimitation. This cutoff
# asks "are these close enough that their reads can be confused", not "are these
# the same species", and BUSCO markers run several points more identical than
# the genome average, so the genome-level number has to sit well below the
# marker-level one it is standing in for.
DEFAULT_ANI_CUTOFF = 0.80
VALID_LINKAGES = ("single", "complete", "average")


class UnionFind:
	"""Disjoint-set with path compression and union by size."""

	def __init__(self):
		self.parent = {}
		self.size = {}

	def add(self, x):
		if x not in self.parent:
			self.parent[x] = x
			self.size[x] = 1

	def find(self, x):
		self.add(x)
		root = x
		while self.parent[root] != root:
			root = self.parent[root]
		# Path compression
		while self.parent[x] != root:
			self.parent[x], x = root, self.parent[x]
		return root

	def union(self, a, b):
		ra, rb = self.find(a), self.find(b)
		if ra == rb:
			return ra
		if self.size[ra] < self.size[rb]:
			ra, rb = rb, ra
		self.parent[rb] = ra
		self.size[ra] += self.size[rb]
		return ra

	def groups(self):
		out = defaultdict(set)
		for x in self.parent:
			out[self.find(x)].add(x)
		return list(out.values())

	def groups_containing(self, x):
		"""Members of x's current group. Used to size a prospective merge."""
		if x not in self.parent:
			return {x}
		root = self.find(x)
		return {y for y in self.parent if self.find(y) == root}


def _detect_scale(values, sample_size=10000):
	"""
	Return the divisor needed to put ANI on a 0-1 scale.

	Distinguishing 0-1 from 0-100 matters: silently treating a 0-100 file as
	0-1 makes every edge pass the cutoff and collapses the entire database
	into one cluster.
	"""
	head = values[:sample_size]
	if not head:
		return 1.0
	if max(head) > 1.0:
		return 100.0
	return 1.0


def load_ani(path, cutoff=DEFAULT_ANI_CUTOFF, restrict_to=None):
	"""
	Read an ANI table and return edges at or above cutoff.

	restrict_to: optional set of genome IDs. Edges touching a genome outside
	this set are dropped, which keeps the graph to genomes the database
	actually contains.

	Returns (edges, stats) where edges is a list of (g1, g2, ani) with ani on
	a 0-1 scale.
	"""
	if not 0.0 < cutoff <= 1.0:
		raise ValueError(
			f"ANI cutoff must be in (0, 1], got {cutoff}. "
			"Cutoffs are expressed on the 0-1 scale (e.g. 0.80)."
		)

	raw = []
	genomes_in_file = set()
	stats = {
		"lines": 0,
		"malformed": 0,
		"self_edges": 0,
		"out_of_db": 0,
		"below_cutoff": 0,
		"kept": 0,
	}

	try:
		handle = open(path)
	except FileNotFoundError:
		raise FileNotFoundError(f"ANI file not found: {path}")

	with handle as f:
		for lineno, line in enumerate(f, 1):
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			stats["lines"] += 1
			parts = line.split("\t")
			if len(parts) < 3:
				stats["malformed"] += 1
				if stats["malformed"] <= 5:
					logger.warning(
						f"{path}:{lineno}: expected 3 tab-separated columns, got "
						f"{len(parts)}; skipping"
					)
				continue
			g1, g2 = parts[0].strip(), parts[1].strip()
			try:
				ani = float(parts[2])
			except ValueError:
				stats["malformed"] += 1
				if stats["malformed"] <= 5:
					logger.warning(f"{path}:{lineno}: non-numeric ANI {parts[2]!r}; skipping")
				continue
			# Membership is recorded before any cutoff or in-database filtering.
			# "Absent from the ANI file" and "present but below cutoff" mean
			# opposite things: the first is an organism that was never measured
			# (transcriptome-derived, no genome to compare), the second is an
			# organism actively determined to have no close relative. Only the
			# first should ever fall back to NCBI taxonomy.
			genomes_in_file.add(g1)
			genomes_in_file.add(g2)
			if g1 == g2:
				stats["self_edges"] += 1
				continue
			raw.append((g1, g2, ani))

	divisor = _detect_scale([a for _, _, a in raw])
	if divisor != 1.0:
		logger.info(f"ANI values look like a 0-100 scale; normalizing to 0-1.")

	edges = []
	for g1, g2, ani in raw:
		ani = ani / divisor
		if restrict_to is not None and (g1 not in restrict_to or g2 not in restrict_to):
			stats["out_of_db"] += 1
			continue
		if ani < cutoff:
			stats["below_cutoff"] += 1
			continue
		edges.append((g1, g2, ani))
		stats["kept"] += 1

	if stats["malformed"] > 5:
		logger.warning(f"{stats['malformed']} malformed lines total in {path}")
	if restrict_to is not None and stats["out_of_db"] and not edges:
		logger.warning(
			"Every ANI edge was dropped as out-of-database. The genome IDs in "
			"the ANI file likely do not match column 3 of the taxid link file."
		)

	# skani (and most ANI tools) only report pairs above an internal screening
	# floor; non-homologous pairs are simply absent rather than listed as low.
	# That is harmless when the floor is well below the cutoff, but if the file
	# was produced with a floor at or above the cutoff, every absent pair is
	# indistinguishable from a genuinely dissimilar one and clustering will
	# silently under-cluster. Surface the observed floor so that is visible.
	if raw:
		observed_floor = min(a / divisor for _, _, a in raw)
		stats["observed_floor"] = observed_floor
		if observed_floor >= cutoff:
			logger.warning(
				f"Lowest ANI in {path} is {observed_floor:.4f}, at or above the "
				f"cutoff of {cutoff}. The file appears to be pre-filtered at or "
				f"above your cutoff, so pairs below it are absent rather than "
				f"low-scoring and cannot be distinguished from unrelated pairs. "
				f"Lowering --ani_cutoff below {observed_floor:.4f} will have no effect."
			)
		elif observed_floor > 0.80:
			logger.info(
				f"Lowest reported ANI is {observed_floor:.4f}; pairs below that "
				f"were not reported by the aligner and are treated as unrelated."
			)

	stats["genomes_in_file"] = genomes_in_file
	logger.info(
		f"ANI: {stats['kept']} edges at >= {cutoff} "
		f"({stats['below_cutoff']} below cutoff, {stats['out_of_db']} out-of-db); "
		f"{len(genomes_in_file)} genomes named in the file"
	)
	return edges, stats


def cluster_genomes(edges, linkage="single", max_component_for_full_linkage=500):
	"""
	Group genomes from ANI edges.

	Returns a list of sets of genome IDs. Genomes with no qualifying edge are
	not returned; callers treat absent genomes as singletons.
	"""
	if linkage not in VALID_LINKAGES:
		raise ValueError(f"linkage must be one of {VALID_LINKAGES}, got {linkage!r}")

	uf = UnionFind()
	for g1, g2, _ in edges:
		uf.union(g1, g2)
	components = uf.groups()

	if linkage == "single":
		return components

	# complete/average linkage refines each single-linkage component.
	sim = {}
	for g1, g2, ani in edges:
		key = (g1, g2) if g1 <= g2 else (g2, g1)
		# Keep the strongest edge if a pair is listed more than once.
		if key not in sim or ani > sim[key]:
			sim[key] = ani

	refined = []
	for comp in components:
		if len(comp) <= 2:
			refined.append(comp)
			continue
		if len(comp) > max_component_for_full_linkage:
			logger.warning(
				f"Component of {len(comp)} genomes exceeds the {linkage}-linkage "
				f"size limit; falling back to single linkage for this component."
			)
			refined.append(comp)
			continue
		refined.extend(_agglomerate(comp, sim, linkage))
	return refined


def _agglomerate(genomes, sim, linkage):
	"""
	Agglomerative clustering within one component. Absent pairs have
	similarity 0, so they can never be merged above a positive cutoff.

	Because this only runs inside an existing single-linkage component and the
	merge criterion is "the two groups are still mutually similar", the result
	is a partition of that component.
	"""
	# Canonical ordering throughout: set iteration order is not stable across
	# interpreter runs, and "first maximum wins" over an unordered list would
	# make cluster assignments vary run to run on identical input.
	groups = [frozenset([g]) for g in sorted(genomes)]

	def pair_sim(a, b):
		vals = []
		for x in a:
			for y in b:
				key = (x, y) if x <= y else (y, x)
				vals.append(sim.get(key, 0.0))
		if not vals:
			return 0.0
		return min(vals) if linkage == "complete" else sum(vals) / len(vals)

	def sort_key(g):
		return sorted(g)

	merged = True
	while merged and len(groups) > 1:
		merged = False
		groups.sort(key=sort_key)
		best = None
		best_score = 0.0
		for i in range(len(groups)):
			for j in range(i + 1, len(groups)):
				s = pair_sim(groups[i], groups[j])
				# Strict > with canonically sorted groups makes the earliest
				# pair win any tie, deterministically.
				if s > best_score:
					best_score, best = s, (i, j)
		# sim already contains only edges >= cutoff, so any positive
		# linkage score means the merge criterion is satisfied.
		if best is not None:
			i, j = best
			groups.append(groups[i] | groups[j])
			for idx in sorted((i, j), reverse=True):
				groups.pop(idx)
			merged = True

	return [set(g) for g in groups]


def load_genome_taxid_map(taxid_link_path):
	"""
	Derive genome <-> species-taxid mappings from busco_taxid_genome_link.txt.

	Columns: marker_seq, taxid, genome_id[,genome_id...]

	Column 3 is comma-separated: one marker can be shared by several genomes.
	The existing eukfrac_calc.py reads this field as a single opaque string,
	which fuses multi-genome fields into one bogus genome name and would break
	any join against an ANI file keyed on real genome IDs.
	"""
	genome_taxids = defaultdict(set)
	taxid_genomes = defaultdict(set)
	seq_genomes = {}
	seq_taxid = {}
	multi_genome_markers = 0
	shared_markers = 0

	with open(taxid_link_path) as f:
		for lineno, line in enumerate(f, 1):
			line = line.rstrip("\n")
			if not line.strip():
				continue
			parts = line.split("\t")
			if len(parts) < 2:
				continue
			seq, taxid = parts[0].strip(), parts[1].strip()
			seq_taxid[seq] = taxid
			if len(parts) < 3 or parts[2].strip() in ("", "NA"):
				seq_genomes[seq] = []
				continue
			genomes = [g.strip() for g in parts[2].split(",") if g.strip() and g.strip() != "NA"]
			seq_genomes[seq] = genomes
			if len(genomes) > 1:
				multi_genome_markers += 1
			if "Collapse" in seq:
				shared_markers += 1

			# Only unambiguous markers establish which species owns which
			# genome. A *Collapse* marker is by construction a sequence that
			# could not be told apart between genomes, and a marker listing
			# several genomes is shared evidence, not ownership. Letting
			# either define ownership chains unrelated species together
			# through shared sequence and silently defeats the ANI cutoff.
			# These markers are still kept in seq_genomes, because reads do
			# align to them and they must be carried into realignment.
			if len(genomes) == 1 and "Collapse" not in seq:
				g = genomes[0]
				genome_taxids[g].add(taxid)
				taxid_genomes[taxid].add(g)

	# A species whose markers are all collapsed/shared would otherwise vanish
	# from the cluster map entirely, so fall back to shared evidence for
	# ownership rather than dropping the taxon.
	for seq, genomes in seq_genomes.items():
		taxid = seq_taxid[seq]
		if taxid not in taxid_genomes and genomes:
			for g in genomes:
				genome_taxids[g].add(taxid)
				taxid_genomes[taxid].add(g)
			logger.debug(
				f"Taxid {taxid} has only shared/collapsed markers; "
				f"using them for genome ownership."
			)

	if multi_genome_markers:
		logger.info(
			f"{multi_genome_markers} markers map to more than one genome "
			f"(comma-separated column 3); treated as shared evidence."
		)
	if shared_markers:
		logger.info(f"{shared_markers} collapsed markers excluded from genome ownership.")
	logger.info(
		f"Loaded {len(seq_taxid)} markers, {len(genome_taxids)} genomes, "
		f"{len(taxid_genomes)} taxa from {taxid_link_path}"
	)
	return {
		"genome_taxids": dict(genome_taxids),
		"taxid_genomes": dict(taxid_genomes),
		"seq_genomes": seq_genomes,
		"seq_taxid": seq_taxid,
	}


def find_uncovered_taxids(taxid_genomes, genomes_in_file):
	"""
	Split taxa by whether the ANI file says anything about them at all.

	This is not the same question as "does this taxon have edges above the
	cutoff". A genome-derived species with no qualifying edge has been measured
	and found to have no close relative, and should stay a singleton. A
	transcriptome-derived organism has no genome to compare and is simply absent
	from the file; leaving it a singleton means it competes with nothing, so if
	it is closely related to a species that really is present it will collect
	that species' reads unchallenged and be reported as a false positive.

	Returns (covered, uncovered, partial) as sorted lists of taxids. `partial`
	means some but not all of a taxon's genomes are in the ANI file, which
	usually indicates an incomplete ANI run rather than a transcriptome entry.
	"""
	covered, uncovered, partial = [], [], []
	for taxid, genomes in taxid_genomes.items():
		seen = [g for g in genomes if g in genomes_in_file]
		if not seen:
			uncovered.append(taxid)
		elif len(seen) < len(genomes):
			partial.append(taxid)
			covered.append(taxid)
		else:
			covered.append(taxid)

	if uncovered:
		logger.info(
			f"{len(uncovered)} of {len(taxid_genomes)} taxa have no entry in the "
			f"ANI file; they can only be grouped by NCBI taxonomy."
		)
	if partial:
		logger.warning(
			f"{len(partial)} taxa have only some of their genomes in the ANI "
			f"file: {','.join(sorted(partial)[:10])}"
			f"{' ...' if len(partial) > 10 else ''}. These are treated as "
			f"covered, but the ANI run may be incomplete."
		)
	return sorted(covered), sorted(uncovered), sorted(partial)


CACHED_RANKS = ("genus", "family", "order", "class", "phylum", "kingdom")
IDENTITY_CACHE_VERSION = 1


def _file_stamp(path):
	try:
		st = os.stat(path)
		return {"size": st.st_size, "mtime": int(st.st_mtime)}
	except OSError:
		return None


def load_identity_cache(db_dir, taxid_link_path=None):
	"""
	Read <database_dir>/identity_cache.json, or return None if it is unusable.

	Unusable means absent, unreadable, written by a different format version, or
	stale with respect to the marker table. Staleness is checked here rather
	than left to the caller so that editing the database cannot silently produce
	results derived from an older taxonomy.
	"""
	if not db_dir:
		return None
	path = os.path.join(db_dir, IDENTITY_CACHE_FILENAME)
	if not os.path.exists(path):
		return None
	try:
		with open(path) as f:
			cache = json.load(f)
	except Exception as e:
		logger.warning(f"Could not read {path} ({e}); it will be rebuilt.")
		return None
	if cache.get("format_version") != IDENTITY_CACHE_VERSION:
		logger.info(
			f"{path} was written by format version "
			f"{cache.get('format_version')}, expected {IDENTITY_CACHE_VERSION}; "
			f"rebuilding."
		)
		return None
	if taxid_link_path:
		stamp = _file_stamp(taxid_link_path)
		if cache.get("source", {}).get("taxid_link") != stamp:
			logger.info(
				f"{path} is stale: {os.path.basename(taxid_link_path)} has "
				f"changed since it was written. Rebuilding."
			)
			return None
	return cache


def build_identity_cache(db_dir, taxid_link_path, ncbi, taxids=None):
	"""
	Build <database_dir>/identity_cache.json and return the payload.

	Resolving a lineage per taxon takes about a minute for a full database and
	the answer depends only on the database, so it is done once and written
	beside the other database files. Called automatically when the cache is
	missing or stale, so a database edit does not require remembering to refresh
	it by hand.

	Rank lineages are cached rather than finished competition groups. Groups
	depend on ani_cutoff, ani_linkage and taxonomy_rank; caching them would go
	stale whenever any of those changed. Rank lineages depend only on the
	taxonomy and the marker table.

	Returns the payload even when it cannot be written, so a read-only database
	directory costs time rather than correctness.
	"""
	from datetime import datetime, timezone

	if taxids is None:
		taxids = set()
		with open(taxid_link_path) as f:
			for line in f:
				parts = line.rstrip("\n").split("\t")
				if len(parts) >= 2 and parts[1].strip():
					taxids.add(parts[1].strip())
		taxids = sorted(taxids)

	logger.info(
		f"Building {IDENTITY_CACHE_FILENAME} for {len(taxids)} taxa. This runs "
		f"once per database and takes about a minute."
	)

	ranks, unresolved, merged = {}, [], {}
	for taxid in taxids:
		try:
			lineage = ncbi.get_lineage(int(taxid))
			if not lineage:
				unresolved.append(taxid)
				continue
			ranks_of = ncbi.get_rank(lineage)
			entry = {}
			for node in lineage:
				r = ranks_of.get(node)
				if r in CACHED_RANKS and r not in entry:
					entry[r] = str(node)
			ranks[taxid] = entry
			if str(lineage[-1]) != taxid:
				merged[taxid] = str(lineage[-1])
		except Exception as e:
			logger.debug(f"no lineage for {taxid}: {e}")
			unresolved.append(taxid)

	payload = {
		"format_version": IDENTITY_CACHE_VERSION,
		"created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
		"cached_ranks": list(CACHED_RANKS),
		"source": {"taxid_link": _file_stamp(taxid_link_path)},
		"n_taxa": len(taxids),
		"ranks": ranks,
		"merged_taxids": merged,
		"unresolved": unresolved,
	}

	path = os.path.join(db_dir, IDENTITY_CACHE_FILENAME) if db_dir else None
	if path and os.access(db_dir, os.W_OK):
		# Unique temp name: samples run in parallel and may build concurrently.
		tmp = f"{path}.tmp.{os.getpid()}"
		try:
			with open(tmp, "w") as f:
				json.dump(payload, f)
			os.replace(tmp, path)
			logger.info(
				f"Wrote {path}: {len(ranks)} taxa, {len(unresolved)} without a "
				f"lineage. Later runs on this database reuse it."
			)
		except Exception as e:
			logger.warning(f"Could not write {path} ({e}); continuing in memory.")
			try:
				os.remove(tmp)
			except OSError:
				pass
	elif path:
		logger.warning(
			f"Database directory is not writable, so {IDENTITY_CACHE_FILENAME} "
			f"cannot be saved and lineages will be resolved again on every run. "
			f"Make {db_dir} writable, or run "
			f"build_db/make_identity_cache.py once as a user who can write there."
		)

	return payload


def rank_map(taxids, ncbi=None, rank="genus", cache=None, db_dir=None,
			 taxid_link_path=None):
	"""
	Map taxids to their NCBI ancestor at `rank`.

	Reads the precomputed cache when one is available, since resolving a lineage
	per taxon takes about a minute for a full database and the answer depends
	only on the database. Falls back to querying the taxonomy for any taxon the
	cache does not cover, so a partially stale cache degrades rather than
	breaking.
	"""
	if cache is None:
		cache = load_identity_cache(db_dir, taxid_link_path)
		if cache is None and ncbi is not None and taxid_link_path:
			# Missing or stale. Build it now so this is a one-off cost rather
			# than something the user has to remember to do.
			cache = build_identity_cache(db_dir, taxid_link_path, ncbi,
										 taxids=list(taxids))

	out = {}
	missing = []

	if cache:
		if rank not in cache.get("cached_ranks", []):
			logger.warning(
				f"{IDENTITY_CACHE_FILENAME} does not contain rank '{rank}' "
				f"(has {','.join(cache.get('cached_ranks', []))}); resolving directly."
			)
		else:
			ranks = cache.get("ranks", {})
			for taxid in taxids:
				entry = ranks.get(taxid)
				if entry is None:
					missing.append(taxid)
				elif rank in entry:
					out[taxid] = entry[rank]
			logger.info(
				f"Rank lineages: {len(out)} taxa from {IDENTITY_CACHE_FILENAME}, "
				f"{len(missing)} not in the cache."
			)
			if not missing:
				return out
	else:
		missing = list(taxids)
		logger.warning(
			f"No usable {IDENTITY_CACHE_FILENAME} and it could not be built; "
			f"resolving lineages directly."
		)

	if not missing:
		return out
	if ncbi is None:
		logger.warning(
			f"{len(missing)} taxa are absent from the cache and no taxonomy was "
			f"supplied; they will not be grouped by {rank}."
		)
		return out

	unranked = []
	for taxid in missing:
		try:
			lineage = ncbi.get_lineage(int(taxid))
			ranks_of = ncbi.get_rank(lineage)
			found = None
			for node in lineage:
				if ranks_of.get(node) == rank:
					found = str(node)
					break
			if found:
				out[taxid] = found
			else:
				unranked.append(taxid)
		except Exception as e:
			logger.debug(f"Could not resolve {rank} for taxid {taxid}: {e}")
			unranked.append(taxid)

	if unranked:
		logger.warning(
			f"{len(unranked)} taxa have no {rank} in their NCBI lineage and stay "
			f"singletons: {','.join(sorted(unranked)[:10])}"
			f"{' ...' if len(unranked) > 10 else ''}"
		)
	return out


# Kept under the old name so external callers do not break.
genus_map_from_ncbi = rank_map


def apply_taxonomy_grouping(uf, taxid_rank_group, restrict_to=None):
	"""
	Union every taxon with its NCBI rank-mates.

	This runs over the whole database, not just taxa missing from the ANI file.
	ANI and taxonomy each catch what the other misses, and the union is more
	inclusive than either:

	  - ANI links near-identical species that NCBI splits across genera, which
	    is the case genus alone cannot see.
	  - Genus links species that ANI cannot measure or reports below the cutoff:
	    transcriptome-derived entries with no genome, pairs below skani's
	    reporting floor, and genuinely divergent congeners that should still be
	    checked against each other.

	Grouping two species together is a decision to *compete* them, not to merge
	them. Divergent members keep their own unique reads and are still called
	separately, so an over-inclusive group costs realignment time rather than
	sensitivity. An under-inclusive one costs a species. That asymmetry is what
	makes unioning the right operation rather than intersecting.

	Mutates `uf` in place. Returns the set of taxids that gained a link here,
	which feeds the provenance label on each cluster.
	"""
	by_group = defaultdict(list)
	for taxid, group_key in taxid_rank_group.items():
		if restrict_to is not None and taxid not in restrict_to:
			continue
		by_group[group_key].append(taxid)

	linked = set()
	for group_key, members in sorted(by_group.items()):
		members = sorted(members)
		if len(members) < 2:
			continue
		for t in members[1:]:
			uf.union(members[0], t)
		linked.update(members)

	multi = sum(1 for m in by_group.values() if len(m) > 1)
	logger.info(
		f"Taxonomy grouping linked {len(linked)} taxa across {multi} groups."
	)
	return linked


def build_species_clusters(genome_clusters, genome_taxids, taxid_genomes,
						   taxid_rank_group=None):
	"""
	Convert genome clusters into species competition groups.

	Three closures are applied, all unioned together:
	  1. genomes -> the species they belong to
	  2. species -> any other cluster holding another genome of that species
	  3. species -> its NCBI rank-mates (when taxid_rank_group is supplied)

	Step 2 matters because a species with several representative genomes could
	otherwise be split across groups and end up competing against itself.
	Step 3 is applied to every taxon, not only those absent from the ANI file;
	see apply_taxonomy_grouping for why union rather than intersection.

	Returns (cluster_id -> set(taxid), taxid -> cluster_id, cluster_id ->
	provenance). Provenance records which evidence linked a group's members:
	"ani", "genus", "ani+genus", or "singleton".
	"""
	uf = UnionFind()
	for taxid in taxid_genomes:
		uf.add(taxid)

	ani_linked = set()
	for comp in genome_clusters:
		taxa = set()
		for g in comp:
			taxa.update(genome_taxids.get(g, ()))
		taxa = sorted(taxa)
		if len(taxa) > 1:
			ani_linked.update(taxa)
		for t in taxa[1:]:
			uf.union(taxa[0], t)

	taxonomy_linked = set()
	if taxid_rank_group:
		taxonomy_linked = apply_taxonomy_grouping(
			uf, taxid_rank_group, restrict_to=set(taxid_genomes)
		)

	clusters = defaultdict(set)
	for taxid in taxid_genomes:
		clusters[uf.find(taxid)].add(taxid)

	# Stable, deterministic cluster IDs so output is reproducible and
	# diffable across runs.
	cluster_map = {}
	taxid_cluster = {}
	provenance = {}
	for i, members in enumerate(sorted(clusters.values(), key=lambda s: sorted(s)), 1):
		cid = f"ANIC{i:06d}"
		cluster_map[cid] = members
		for t in members:
			taxid_cluster[t] = cid
		if len(members) == 1:
			provenance[cid] = "singleton"
			continue
		has_ani = bool(members & ani_linked)
		has_tax = bool(members & taxonomy_linked)
		if has_ani and has_tax:
			provenance[cid] = "ani+genus"
		elif has_ani:
			provenance[cid] = "ani"
		else:
			provenance[cid] = "genus"

	sizes = sorted((len(m) for m in cluster_map.values()), reverse=True)
	multi = sum(1 for n in sizes if n > 1)
	by_prov = defaultdict(int)
	for p in provenance.values():
		by_prov[p] += 1
	logger.info(
		f"Built {len(cluster_map)} competition groups ({multi} with more than "
		f"one species): " + ", ".join(f"{k}={v}" for k, v in sorted(by_prov.items()))
	)
	if sizes and sizes[0] > 1:
		logger.info(
			f"Largest groups: {sizes[:5]}. Realignment cost depends on how many "
			f"members have reads in a given sample, not on group size."
		)
	return cluster_map, taxid_cluster, provenance


def build(ani_path, taxid_link_path, cutoff=DEFAULT_ANI_CUTOFF, linkage="single",
		  ncbi=None, taxonomy_rank="genus", group_by_taxonomy=True,
		  db_dir=None):
	"""
	ANI table + taxid link file -> species competition groups.

	Groups are the union of two kinds of link: genomes at or above the ANI
	cutoff, and taxa sharing an NCBI ancestor at `taxonomy_rank`. Neither alone
	is sufficient. ANI cannot see transcriptome-derived entries or pairs below
	the aligner's reporting floor; NCBI genus cannot see near-identical species
	that taxonomy has split across genera. Taking the union means a wrong or
	over-split genus is corrected by ANI, and a gap in ANI coverage is covered
	by taxonomy.

	Rank lineages come from <db_dir>/identity_cache.json when it exists.
	Resolving them from the taxonomy takes about a minute for a full database
	and would otherwise repeat once per sample. Only the lineages are cached,
	not the finished groups: groups depend on the cutoff, the linkage and the
	rank, so caching them would go stale whenever any of those changed, while
	the clustering itself takes milliseconds.
	"""
	if db_dir is None:
		db_dir = os.path.dirname(os.path.abspath(taxid_link_path))

	link = load_genome_taxid_map(taxid_link_path)
	edges, stats = load_ani(
		ani_path, cutoff=cutoff, restrict_to=set(link["genome_taxids"])
	)
	genome_clusters = cluster_genomes(edges, linkage=linkage)

	covered, uncovered, partial = find_uncovered_taxids(
		link["taxid_genomes"], stats.get("genomes_in_file", set())
	)

	rank_group = None
	if group_by_taxonomy:
		rank_group = rank_map(
			list(link["taxid_genomes"]), ncbi=ncbi, rank=taxonomy_rank,
			db_dir=db_dir, taxid_link_path=taxid_link_path,
		)
		if not rank_group:
			logger.warning(
				f"No {taxonomy_rank} could be resolved for any taxon, so groups "
				f"come from ANI alone."
			)

	cluster_map, taxid_cluster, provenance = build_species_clusters(
		genome_clusters, link["genome_taxids"], link["taxid_genomes"],
		taxid_rank_group=rank_group,
	)

	if uncovered:
		still_alone = [t for t in uncovered if len(cluster_map[taxid_cluster[t]]) == 1]
		logger.info(
			f"{len(uncovered)} taxa have no ANI entry; "
			f"{len(uncovered) - len(still_alone)} were grouped by taxonomy, "
			f"{len(still_alone)} remain unable to compete with anything."
		)

	return {
		"clusters": cluster_map,
		"taxid_cluster": taxid_cluster,
		"provenance": provenance,
		"link": link,
		"ani_stats": {k: v for k, v in stats.items() if k != "genomes_in_file"},
		"covered": covered,
		"uncovered": uncovered,
		"partial": partial,
	}


def main(argv=None):
	parser = argparse.ArgumentParser(
		description="Build ANI-based species competition groups for EukDetect."
	)
	parser.add_argument("--ani", required=True, help="genome1<TAB>genome2<TAB>ANI")
	parser.add_argument("--taxid_link", required=True, help="busco_taxid_genome_link.txt")
	parser.add_argument(
		"--ani_cutoff", type=float, default=DEFAULT_ANI_CUTOFF,
		help=f"ANI cutoff on the 0-1 scale (default {DEFAULT_ANI_CUTOFF})",
	)
	parser.add_argument(
		"--linkage", choices=VALID_LINKAGES, default="single",
		help="Linkage criterion (default single)",
	)
	parser.add_argument(
		"--taxdb", default=None,
		help="ete3 NCBI taxonomy sqlite. Required for taxonomy grouping.",
	)
	parser.add_argument(
		"--taxonomy_rank", default="genus",
		help="NCBI rank unioned with ANI to form groups (default genus)",
	)
	parser.add_argument(
		"--no_taxonomy_grouping", action="store_true",
		help="Use ANI links only, without unioning in NCBI rank-mates",
	)
	parser.add_argument(
		"--database_dir", default=None,
		help="Database directory holding identity_cache.json. Defaults to the "
			 "directory containing --taxid_link.",
	)
	parser.add_argument("--out", default="-", help="Output TSV, or - for stdout")
	args = parser.parse_args(argv)

	logging.basicConfig(
		format="%(asctime)s [%(levelname)s] %(message)s",
		datefmt="%Y-%m-%d %H:%M:%S",
		level=logging.INFO,
	)

	ncbi = None
	if args.taxdb and not args.no_taxonomy_grouping:
		from ete3 import NCBITaxa
		ncbi = NCBITaxa(args.taxdb)

	result = build(
		args.ani, args.taxid_link, args.ani_cutoff, args.linkage,
		ncbi=ncbi, taxonomy_rank=args.taxonomy_rank,
		group_by_taxonomy=not args.no_taxonomy_grouping,
		db_dir=args.database_dir,
	)

	dest = sys.stdout if args.out == "-" else open(args.out, "w")
	try:
		dest.write("cluster_id\tn_species\tprovenance\ttaxids\tgenomes\n")
		uncovered = set(result["uncovered"])
		for cid, members in sorted(result["clusters"].items()):
			genomes = set()
			for t in members:
				genomes.update(result["link"]["taxid_genomes"].get(t, ()))
			prov = result["provenance"].get(cid, "ani")
			if len(members) == 1 and members & uncovered:
				prov = "uncovered_singleton"
			dest.write(
				f"{cid}\t{len(members)}\t{prov}\t"
				f"{','.join(sorted(members))}\t{','.join(sorted(genomes))}\n"
			)
	finally:
		if dest is not sys.stdout:
			dest.close()


if __name__ == "__main__":
	main()
