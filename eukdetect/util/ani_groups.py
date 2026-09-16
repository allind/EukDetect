#!/usr/bin/env python

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)


DEFAULT_ANI_THRESHOLD = 80.0
class _Union:
	"""Disjoint-set over taxids."""

	def __init__(self):
		self.parent = {}

	def add(self, x):
		self.parent.setdefault(x, x)

	def find(self, x):
		self.add(x)
		root = x
		while self.parent[root] != root:
			root = self.parent[root]
		while self.parent[x] != root:
			self.parent[x], x = root, self.parent[x]
		return root

	def union(self, a, b):
		ra, rb = self.find(a), self.find(b)
		if ra != rb:
			self.parent[rb] = ra


def load_ani_pairs(path, threshold=DEFAULT_ANI_THRESHOLD):

	pairs = []
	stats = {"lines": 0, "malformed": 0, "self": 0, "below": 0, "kept": 0,
			 "min_reported": None, "scaled": False}
	raw = []

	try:
		with open(path) as f:
			for lineno, line in enumerate(f, 1):
				line = line.strip()
				if not line or line.startswith("#"):
					continue
				stats["lines"] += 1
				parts = line.split("\t")
				if len(parts) < 3:
					stats["malformed"] += 1
					continue
				g1, g2 = parts[0].strip(), parts[1].strip()
				try:
					ani = float(parts[2])
				except ValueError:
					stats["malformed"] += 1
					continue
				if g1 == g2:
					stats["self"] += 1
					continue
				raw.append((g1, g2, ani))
	except FileNotFoundError:
		logger.error(f"ANI file not found: {path}")
		return [], stats
	except Exception as e:
		logger.error(f"Could not read ANI file {path}: {e}")
		return [], stats

	if not raw:
		logger.warning(f"No usable pairs in {path}")
		return [], stats

	# A 0-1 file compared against a 0-100 threshold would drop every pair, and a
	# 0-100 file against a 0-1 threshold would keep every pair. Detect which.
	if max(a for _, _, a in raw) <= 1.0:
		raw = [(g1, g2, a * 100.0) for g1, g2, a in raw]
		stats["scaled"] = True
		logger.info("ANI values look like a 0-1 scale; read as percentages.")

	stats["min_reported"] = min(a for _, _, a in raw)
	for g1, g2, ani in raw:
		if ani >= threshold:
			pairs.append((g1, g2))
			stats["kept"] += 1
		else:
			stats["below"] += 1

	if stats["malformed"]:
		logger.warning(f"{stats['malformed']} malformed lines in {path}")
	logger.info(
		f"ANI: {stats['kept']} genome pairs at or above {threshold} "
		f"({stats['below']} below) from {stats['lines']} lines"
	)
	# Aligners only report pairs above an internal floor. If that floor sits
	# above the threshold the threshold is doing nothing, which is worth saying
	# plainly rather than leaving it to look like a setting that had an effect.
	if stats["min_reported"] is not None and stats["min_reported"] >= threshold:
		logger.info(
			f"Lowest ANI reported is {stats['min_reported']:.2f}, at or above "
			f"the {threshold} threshold, so every reported pair is being used."
		)
	return pairs, stats


def group_species(genuses, genome_taxids, ani_pairs):

	uf = _Union()
	genus_of = {}

	for genus, taxids in genuses.items():
		for t in taxids:
			uf.add(t)
			genus_of[t] = genus
		for t in taxids[1:]:
			uf.union(taxids[0], t)

	linked_by_ani = 0
	cross_genus = 0
	for g1, g2 in ani_pairs:
		t1 = genome_taxids.get(g1)
		t2 = genome_taxids.get(g2)
		if not t1 or not t2:
			continue
		for a in t1:
			for b in t2:
				if a == b:
					continue
				# Only taxa the sample actually contains are in genus_of; ANI
				# covers the whole database, most of which is irrelevant here.
				if a not in genus_of or b not in genus_of:
					continue
				if uf.find(a) != uf.find(b):
					linked_by_ani += 1
					if genus_of.get(a) != genus_of.get(b):
						cross_genus += 1
					uf.union(a, b)

	grouped = defaultdict(list)
	for t in genus_of:
		grouped[uf.find(t)].append(t)

	groups = {}
	for i, (_root, members) in enumerate(sorted(grouped.items(),
											   key=lambda kv: sorted(kv[1])), 1):
		groups[f"GRP{i:05d}"] = sorted(members)

	stats = {
		"n_groups": len(groups),
		"n_multi": sum(1 for m in groups.values() if len(m) > 1),
		"ani_merges": linked_by_ani,
		"cross_genus_merges": cross_genus,
	}
	logger.info(
		f"Disambiguation groups: {len(groups)} "
		f"({stats['n_multi']} with more than one species). "
		f"ANI merged {linked_by_ani} pairs that genus alone would have "
		f"separated, {cross_genus} of them across different genera."
	)
	return groups, stats
