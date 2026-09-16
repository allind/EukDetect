#!/usr/bin/env python
import logging
import math
import random
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

DEFAULT_MIN_READS = 4
DEFAULT_MIN_MARKERS = 2

DEFAULT_ALPHA = 0.05


DEFAULT_SMALL_N_THRESHOLD = 15


DEFAULT_SMALL_N_MAX_GAP = 1.5


DEFAULT_BASE_MIN_GAP = 0.75
DEFAULT_ABUNDANCE_SCALE = 1.5

DEFAULT_ANCHOR_MARGIN = 1.5

DEFAULT_SPLIT_N_BOOT = 500
DEFAULT_SPLIT_ALPHA = 0.01


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


def _mad(values):
	"""Median absolute deviation: spread unmoved by a handful of outliers."""
	med = _median(values)
	if med is None:
		return 0.0
	return _median([abs(v - med) for v in values]) or 0.0


_MIN_POOL_RATIO = 1


def split_gap(values):

	v = sorted(values)
	n = len(v)
	if n < 4:
		return 0.0
	mid = n // 2
	return _median(v[n - mid:]) - _median(v[:mid])


def bootstrap_split_pvalue(sec, pri, n_boot=1000, rng=None):

	sec_ids = [i for i, _ in sec]
	obs = split_gap(sec_ids)
	n = len(sec_ids)
	if n < 4 or len(pri) < n * _MIN_POOL_RATIO:
		return None, obs

	rng = rng or random
	pri_by_busco = defaultdict(list)
	for ident, marker in pri:
		pri_by_busco[marker_busco(marker)].append(ident)
	pri_pool_all = [i for i, _ in pri]
	sec_buscos = [marker_busco(m) for _, m in sec]

	ge = 0
	for _ in range(n_boot):
		sample = []
		for busco in sec_buscos:

			pool = pri_by_busco.get(busco) or pri_pool_all
			sample.append(rng.choice(pool))
		if split_gap(sample) >= obs:
			ge += 1

	return (ge + 1) / (n_boot + 1), obs


def mwu_less_pvalue(x, y):

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
			 abundance_scale=DEFAULT_ABUNDANCE_SCALE,
			 anchor_margin=DEFAULT_ANCHOR_MARGIN,
			 anchor_alpha=DEFAULT_SPLIT_ALPHA):

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
		"anchor_reads": 0,
		"anchor_markers": 0,
		"anchor_split_gap": None,
		"anchor_p_value": None,
		"anchor_override": False,
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
	pri = []  # (identity, marker) pairs, kept for the stratified bootstrap below
	for p in primaries:
		pri.extend(by_taxid.get(p, []))
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

	ratio = len(sec_ids) / len(pri_ids) if pri_ids else 0.0
	required_gap = base_min_gap + abundance_scale * min(ratio, 1.0)
	rec["abundance_ratio"] = ratio
	rec["required_gap"] = required_gap

	if p is not None and p < alpha and gap >= required_gap:

		anchor_p, split_obs = bootstrap_split_pvalue(sec, pri)
		anchor_threshold = pri_med - anchor_margin
		anchor = [(i, m) for i, m in sec if i >= anchor_threshold]
		anchor_buscos = {marker_busco(m) for _, m in anchor}
		rec["anchor_reads"] = len(anchor)
		rec["anchor_markers"] = len(anchor_buscos)
		rec["anchor_split_gap"] = split_obs
		rec["anchor_p_value"] = anchor_p

		enough_excess = anchor_p is not None and anchor_p < anchor_alpha

		enough_evidence = len(anchor) >= min_reads and len(anchor_buscos) >= min_markers

		if enough_excess and enough_evidence:
			rec["promote"] = True
			rec["anchor_override"] = True
			rec["reason"] = (
				f"kept: reads match it at {sec_med:.2f}% overall (which alone "
				f"would be reassigned, {gap:.2f} points below the primary's "
				f"{pri_med:.2f}%, Mann-Whitney p={p:.2e}), but its own reads split "
				f"into a much wider high/low gap ({split_obs:.2f} points) than "
				f"same-sized draws from the primary typically show (bootstrap "
				f"p={anchor_p:.3f}), with {len(anchor)} reads across "
				f"{len(anchor_buscos)} markers reaching within {anchor_margin} "
				f"points of the primary's level -- consistent with a real presence "
				f"diluted by spillover rather than pure spillover alone"
			)
			return rec

		rec["reason"] = (
			f"reassigned: reads match it at {sec_med:.2f}%, significantly "
			f"below the primary's {pri_med:.2f}% by {gap:.2f} points, at or "
			f"above the {required_gap:.2f}-point minimum required given it is "
			f"{ratio:.1%} as abundant as its primary (Mann-Whitney p={p:.2e}, "
			f"n={len(sec_ids)} vs {len(pri_ids)}; its own split-gap of "
			f"{split_obs:.2f} points did not look bigger than ordinary noise "
			f"produces for a real single population of this size, bootstrap "
			f"p={'NA' if anchor_p is None else format(anchor_p, '.3f')})"
		)
		return rec
	if p is not None and p < alpha:

		logger.debug(
			f"{secondary}: significant (p={p:.2e}) but gap {gap:.2f} is below "
			f"the {required_gap:.2f}-point minimum required at {ratio:.1%} "
			f"relative abundance; not reassigning on significance alone."
		)


	if len(sec_ids) < small_n_threshold and gap > small_n_max_gap:
		rec["reason"] = (
			f"reassigned: only {len(sec_ids)} reads, too few for the "
			f"significance test to rule out a real difference "
			f"(p={p:.2f}), and its {gap:.2f}-point gap from the primary is "
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
			"Abundance_ratio\tRequired_gap\t"
			"Anchor_override\tSplit_gap\tSplit_p_value\t"
			"Anchor_reads\tAnchor_markers\tReason\n"
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
				f"{'yes' if rec['anchor_override'] else 'no'}\t"
				f"{fmt(rec['anchor_split_gap'])}\t{fmt(rec['anchor_p_value'], '.3e')}\t"
				f"{rec['anchor_reads']}\t{rec['anchor_markers']}\t"
				f"{rec['reason']}\n"
			)

	kept = sum(1 for r in records if r["promote"])
	logger.info(
		f"Reassignment: {len(records)} secondaries evaluated, {kept} kept on "
		f"their own read evidence, {len(records) - kept} reassigned. "
		f"Details in {path}"
	)
