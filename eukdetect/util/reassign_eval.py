#!/usr/bin/env python
"""Read-level re-evaluation of proposed species reassignments.

Earlier stages of eukfrac_calc propose that a species is spillover from a
competitor ("primary") and should have its reads folded into it. This module
re-tests every such proposal against the per-read identities in the filtered
BAM, and overturns the proposal when the evidence does not support it.

The decision is two parts, both computed from per-read percent identity:

  1. A one-sided Mann-Whitney U test asking whether the secondary's reads sit
     stochastically lower than its primary's.
  2. A minimum difference between the two medians, scaled by how abundant the
     secondary is relative to its primary, measured in reads per base of marker
     rather than raw reads.

A secondary is reassigned only when both fire. A secondary with too few reads
for the test to have power is reassigned if its median sits far enough below.
Anything else is kept as a separate species.
"""
import logging
import math
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

#A secondary needs at least this much evidence to be judged on its reads at all.
#Below either threshold the proposed reassignment stands unexamined.
DEFAULT_MIN_READS = 4
DEFAULT_MIN_MARKERS = 2

DEFAULT_ALPHA = 0.05

#Below this many reads the Mann-Whitney test has little power, so a large
#median difference is taken as sufficient evidence on its own.
DEFAULT_SMALL_N_THRESHOLD = 15
DEFAULT_SMALL_N_MAX_GAP = 1.5

#required_gap = base + scale * min(abundance_ratio, 1), in percentage points of
#identity. The floor applies to a vanishingly rare secondary; the ceiling
#(base + scale) to one as dense as its primary. Scaling with abundance is what
#keeps genuinely co-present congeners apart: two species at similar depth need a
#much larger identity difference before one is called spillover from the other.
DEFAULT_BASE_MIN_GAP = 0.75
DEFAULT_ABUNDANCE_SCALE = 1.5


BUSCO_RE = re.compile(r"-(\d+at\d+)-")


def marker_busco(marker):

	m = BUSCO_RE.search(marker)
	return m.group(1) if m else marker


def _median(values):
	if not values:
		return None
	v = sorted(values)
	n = len(v)
	return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def _abundance_ratio(secondary, primaries, n_sec, n_pri, marker_lengths):
	# Reads per base of marker, not raw reads: a database genome that carries
	# fewer markers collects fewer reads at the same coverage, which would
	# otherwise make a real species look scarcer than it is.
	if marker_lengths:
		sec_len = marker_lengths.get(secondary, 0)
		pri_len = sum(marker_lengths.get(p, 0) for p in primaries)
		if sec_len > 0 and pri_len > 0 and n_pri > 0:
			return (n_sec / sec_len) / (n_pri / pri_len)
	return n_sec / n_pri if n_pri else 0.0


def mwu_less_pvalue(x, y):
	"""One-sided Mann-Whitney U, normal approximation with tie and continuity
	correction. Returns P(x stochastically below y), or None when undefined
	(an empty group, or every value tied)."""
	n1, n2 = len(x), len(y)
	if n1 == 0 or n2 == 0:
		return None

	combined = [(v, 0) for v in x] + [(v, 1) for v in y]
	combined.sort(key=lambda t: t[0])

	ranks = [0.0] * len(combined)
	tie_term = 0
	i = 0
	while i < len(combined):
		j = i
		while j < len(combined) and combined[j][0] == combined[i][0]:
			j += 1
		avg_rank = (i + 1 + j) / 2.0
		for k in range(i, j):
			ranks[k] = avg_rank
		t = j - i
		if t > 1:
			tie_term += t ** 3 - t
		i = j

	R1 = sum(r for (_, g), r in zip(combined, ranks) if g == 0)
	U1 = R1 - n1 * (n1 + 1) / 2.0
	N = n1 + n2
	if N <= 1:
		return None

	mean_U = n1 * n2 / 2.0
	var_U = n1 * n2 * (N + 1) / 12.0 - n1 * n2 * tie_term / (12.0 * N * (N - 1))
	if var_U <= 0:
		return None
	sd = math.sqrt(var_U)

	cc = 0.5 if U1 < mean_U else (-0.5 if U1 > mean_U else 0.0)
	z = (U1 - mean_U + cc) / sd
	return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def describe_parameters(alpha=DEFAULT_ALPHA, base_min_gap=DEFAULT_BASE_MIN_GAP,
						abundance_scale=DEFAULT_ABUNDANCE_SCALE,
						min_reads=DEFAULT_MIN_READS, min_markers=DEFAULT_MIN_MARKERS):
	"""One line for the log recording the thresholds this run used, so a result
	can be tied to its settings without reading the source."""
	return (
		f"Reassignment test: Mann-Whitney alpha={alpha:g}, median difference must "
		f"reach {base_min_gap:g} + {abundance_scale:g} x min(abundance ratio, 1) "
		f"points ({base_min_gap:g} for a vanishingly rare secondary, "
		f"{base_min_gap + abundance_scale:g} at parity); abundance ratio is reads "
		f"per base of marker; a secondary needs >={min_reads} reads on "
		f">={min_markers} markers to be judged on its reads at all."
	)


def load_read_evidence(path, seq_taxid):

	by_taxid = defaultdict(list)
	skipped_unknown_taxon = 0
	skipped_malformed = 0

	try:
		with open(path) as f:
			header = f.readline().rstrip("\n").split("\t")
			idx = {name: i for i, name in enumerate(header)}
			required = ("Marker", "Identity")
			missing = [c for c in required if c not in idx]
			if missing:
				logger.error(
					f"{path} is missing column(s) {missing}; found {header}. "
					f"Reassignment cannot proceed without per-read identity."
				)
				return {}

			for line in f:
				r = line.rstrip("\n").split("\t")
				if len(r) < len(header):
					skipped_malformed += 1
					continue
				marker = r[idx["Marker"]]
				taxid = seq_taxid.get(marker)
				if taxid is None:
					skipped_unknown_taxon += 1
					continue
				try:
					ident = float(r[idx["Identity"]])
				except ValueError:
					skipped_malformed += 1
					continue
				by_taxid[taxid].append((ident, marker))
	except FileNotFoundError:
		logger.error(f"Per-read evidence not found: {path}")
		return {}

	if skipped_unknown_taxon:
		logger.debug(
			f"{skipped_unknown_taxon} read alignments were on markers not in "
			f"the taxid link table and were skipped."
		)
	if skipped_malformed:
		logger.warning(f"{skipped_malformed} malformed rows in {path} were skipped.")
	logger.info(
		f"Loaded per-read evidence for {len(by_taxid)} taxa, "
		f"{sum(len(v) for v in by_taxid.values())} reads, from {path}"
	)
	return dict(by_taxid)


def evaluate(secondary, primaries, by_taxid,
			 min_reads=DEFAULT_MIN_READS, min_markers=DEFAULT_MIN_MARKERS,
			 alpha=DEFAULT_ALPHA, small_n_threshold=DEFAULT_SMALL_N_THRESHOLD,
			 small_n_max_gap=DEFAULT_SMALL_N_MAX_GAP,
			 base_min_gap=DEFAULT_BASE_MIN_GAP,
			 abundance_scale=DEFAULT_ABUNDANCE_SCALE, marker_lengths=None):
	"""Decide whether one proposed reassignment should stand.

	`promote` True means the proposal is overturned and the secondary is kept as
	a separate species; False means the reassignment stands."""
	sec = by_taxid.get(secondary, [])
	sec_ids = [i for i, _ in sec]
	sec_buscos = {marker_busco(m) for _, m in sec}

	rec = {
		"secondary": secondary,
		"primaries": sorted(primaries),
		"secondary_reads": len(sec_ids),
		"secondary_markers": len(sec_buscos),
		"secondary_median_identity": _median(sec_ids),
		"primary_reads": 0,
		"primary_median_identity": None,
		"p_value": None,
		"gap": None,
		"abundance_ratio": None,
		"required_gap": None,
		"promote": False,
		"reason": "",
	}

	if len(sec_ids) < min_reads:
		rec["reason"] = f"reassigned: only {len(sec_ids)} reads, below the {min_reads} required"
		return rec
	if len(sec_buscos) < min_markers:
		rec["reason"] = (f"reassigned: only {len(sec_buscos)} marker(s), below "
						 f"the {min_markers} required")
		return rec

	pri_ids = []
	for p in primaries:
		pri_ids.extend(i for i, _ in by_taxid.get(p, []))
	rec["primary_reads"] = len(pri_ids)

	if not pri_ids:
		rec["reason"] = "reassigned: no primary reads to compare against"
		return rec

	pri_med = _median(pri_ids)
	sec_med = rec["secondary_median_identity"]
	gap = pri_med - sec_med
	rec["primary_median_identity"] = pri_med
	rec["gap"] = gap

	p = mwu_less_pvalue(sec_ids, pri_ids)
	rec["p_value"] = p

	ratio = _abundance_ratio(secondary, primaries, len(sec_ids), len(pri_ids), marker_lengths)
	required_gap = base_min_gap + abundance_scale * min(ratio, 1.0)
	rec["abundance_ratio"] = ratio
	rec["required_gap"] = required_gap

	#Both conditions must hold: the difference is detectable AND large enough to
	#mean something at this relative abundance.
	if p is not None and p < alpha and gap >= required_gap:
		rec["reason"] = (
			f"reassigned: reads match it at {sec_med:.2f}%, significantly "
			f"below the primary's {pri_med:.2f}% by {gap:.2f} points, at or "
			f"above the {required_gap:.2f}-point minimum required given it is "
			f"{ratio:.1%} as abundant as its primary (Mann-Whitney p={p:.2e}, "
			f"n={len(sec_ids)} vs {len(pri_ids)})"
		)
		return rec

	if p is not None and p < alpha:
		logger.debug(
			f"{secondary}: significant (p={p:.2e}) but gap {gap:.2f} is below "
			f"the {required_gap:.2f}-point minimum required at {ratio:.1%} "
			f"relative abundance; not reassigning on significance alone."
		)

	#Too few reads for the test to rule out a real difference, but the medians
	#are far enough apart that spillover is the better explanation.
	if len(sec_ids) < small_n_threshold and gap > small_n_max_gap:
		rec["reason"] = (
			f"reassigned: only {len(sec_ids)} reads, too few for the "
			f"significance test to rule out a real difference "
			f"(p={'NA' if p is None else format(p, '.2f')}), and its "
			f"{gap:.2f}-point gap from the primary is "
			f"larger than the {small_n_max_gap} points tolerated at this "
			f"sample size"
		)
		return rec

	rec["promote"] = True
	if p is None:
		rec["reason"] = f"kept: reads match it at {sec_med:.2f}%, no test statistic available"
	elif p < alpha:
		rec["reason"] = (
			f"kept: reads match it at {sec_med:.2f}%, {gap:.2f} points below "
			f"the primary's {pri_med:.2f}% -- statistically detectable "
			f"(Mann-Whitney p={p:.2e}) but below the {required_gap:.2f}-point "
			f"minimum required at {ratio:.1%} relative abundance, so not "
			f"treated as a meaningfully different species"
		)
	else:
		rec["reason"] = (
			f"kept: reads match it at {sec_med:.2f}%, not significantly below "
			f"the primary's {pri_med:.2f}% (Mann-Whitney p={p:.2f}, "
			f"n={len(sec_ids)} vs {len(pri_ids)})"
		)
	return rec


def write_report(path, records):

	with open(path, "w") as f:
		f.write(
			"Secondary_taxid\tPrimary_taxids\tPromoted\t"
			"Secondary_reads\tSecondary_markers\tSecondary_median_identity\t"
			"Primary_reads\tPrimary_median_identity\tGap\tMWU_p_value\t"
			"Abundance_ratio\tRequired_gap\tReason\n"
		)
		for rec in records:
			def fmt(v, spec=".3f"):
				return "NA" if v is None else format(v, spec)
			f.write(
				f"{rec['secondary']}\t{','.join(rec['primaries'])}\t"
				f"{'yes' if rec['promote'] else 'no'}\t"
				f"{rec['secondary_reads']}\t{rec['secondary_markers']}\t"
				f"{fmt(rec['secondary_median_identity'])}\t"
				f"{rec['primary_reads']}\t{fmt(rec['primary_median_identity'])}\t"
				f"{fmt(rec['gap'])}\t{fmt(rec['p_value'], '.3e')}\t"
				f"{fmt(rec['abundance_ratio'], '.3f')}\t{fmt(rec['required_gap'])}\t"
				f"{rec['reason']}\n"
			)

	kept = sum(1 for r in records if r["promote"])
	logger.info(
		f"Reassignment: {len(records)} secondaries evaluated, {kept} kept on "
		f"their own read evidence, {len(records) - kept} reassigned. "
		f"Details in {path}"
	)
