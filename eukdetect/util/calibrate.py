#!/usr/bin/env python


import logging
import os
from collections import defaultdict

logger = logging.getLogger(__name__)


DEFAULT_ANCHOR_IDENTITY = 0.97


DEFAULT_ANCHOR_QUANTILE = 0.90


DEFAULT_ANCHOR_SLACK = 0.02

DEFAULT_MIN_SEPARATION = 0.01


DEFAULT_MIN_ANCHOR_IDENTITY = 0.95


ANCHOR_QUANTILE = 0.75

ANCHOR_TOLERANCE = 0.015

DEFAULT_MIN_ANCHOR_READS = 4
DEFAULT_MIN_BRIDGE_READS = 3


def _median(values):
	vals = sorted(values)
	n = len(vals)
	if n == 0:
		return None
	mid = n // 2
	return vals[mid] if n % 2 else (vals[mid - 1] + vals[mid]) / 2.0


def _quantile(values, q):
	if not values:
		return None
	v = sorted(values)
	return v[min(int(q * len(v)), len(v) - 1)]


def estimate_anchor_threshold(fragment_results,
							  anchor_identity=DEFAULT_ANCHOR_IDENTITY,
							  min_anchor_identity=DEFAULT_MIN_ANCHOR_IDENTITY,
							  min_anchor_reads=DEFAULT_MIN_ANCHOR_READS,
							  cluster_id="?"):

	by_best = defaultdict(list)
	for res in fragment_results.values():
		if len(res["winners"]) == 1:
			t = res["winners"][0]
			i = res.get("identity", {}).get(t)
			if i is not None:
				by_best[t].append(i)

	levels = {
		t: _quantile(v, ANCHOR_QUANTILE)
		for t, v in by_best.items()
		if len(v) >= min_anchor_reads
	}
	if not levels:
		return anchor_identity, {}

	observed = max(levels.values())
	threshold = observed - ANCHOR_TOLERANCE
	clamped = min(max(threshold, min_anchor_identity), anchor_identity)

	if clamped != threshold:
		why = ("clamped up to the floor" if clamped > threshold
			   else "capped at the configured maximum")
		logger.info(
			f"{cluster_id}: within-species level estimated at "
			f"{observed:.1%}, giving a threshold of {threshold:.1%}, {why} "
			f"({clamped:.1%})."
		)
	else:
		logger.info(
			f"{cluster_id}: within-species level estimated at {observed:.1%}; "
			f"fragments at or above {clamped:.1%} count as evidence."
		)
	return clamped, levels


def estimate_anchor_threshold(fragment_results,
							  quantile=DEFAULT_ANCHOR_QUANTILE,
							  slack=DEFAULT_ANCHOR_SLACK,
							  floor=DEFAULT_ANCHOR_IDENTITY):

	idents = []
	for res in fragment_results.values():
		if len(res["winners"]) == 1:
			i = res.get("identity", {}).get(res["winners"][0])
			if i is not None:
				idents.append(i)
	if not idents:
		return floor, None, 0

	idents.sort()
	idx = min(int(quantile * len(idents)), len(idents) - 1)
	level = idents[idx]
	return max(level - slack, floor), level, len(idents)


def find_anchors(fragment_results, anchor_identity=DEFAULT_ANCHOR_IDENTITY,
				 min_anchor_reads=DEFAULT_MIN_ANCHOR_READS):

	best_hits = defaultdict(list)
	for res in fragment_results.values():
		winners = res["winners"]
		if len(winners) != 1:
			continue
		t = winners[0]
		ident = res.get("identity", {}).get(t)
		if ident is not None:
			best_hits[t].append(ident)

	anchors = {}
	for t, idents in best_hits.items():
		clean = [i for i in idents if i >= anchor_identity]
		if len(clean) < min_anchor_reads:
			continue
		anchors[t] = _median(clean)
	return anchors


def estimate_delta(fragment_results, anchor, candidate,
				   min_bridge_reads=DEFAULT_MIN_BRIDGE_READS,
				   anchor_identity=DEFAULT_ANCHOR_IDENTITY):

	drops = []
	for res in fragment_results.values():
		idents = res.get("identity", {})
		if anchor in idents and candidate in idents:
			if idents[anchor] < anchor_identity:
				continue
			drops.append(idents[anchor] - idents[candidate])
	if len(drops) < min_bridge_reads:
		return None, len(drops)
	return _median(drops), len(drops)


def classify_fragments(fragment_results, anchors, deltas,
					   min_separation=DEFAULT_MIN_SEPARATION, merged=None,
					   anchor_identity=DEFAULT_ANCHOR_IDENTITY):

	merged = merged or {}

	dominant = None
	if not anchors:
		support = defaultdict(int)
		for res in fragment_results.values():
			if len(res["winners"]) == 1:
				support[res["winners"][0]] += 1
		if support:
			dominant = sorted(support, key=lambda t: (-support[t], t))[0]

	out = {}
	for qname, res in fragment_results.items():
		winners = res["winners"]
		idents = res.get("identity", {})
		if not winners:
			continue
		if len(winners) > 1:
			out[qname] = ("tied", None, None)
			continue
		best = winners[0]

		if best in merged:
			out[qname] = ("undecidable", None, best)
			continue

		obs = idents.get(best)
		if obs is None:
			out[qname] = ("unanchored", best, best)
			continue

		if obs >= anchor_identity:
			out[qname] = ("anchor" if best in anchors else "resolved", best, best)
			continue

		if not anchors:

			target = dominant if (dominant and dominant in idents) else best
			out[qname] = ("unanchored", target, best)
			continue

		target, score = None, None
		if best in anchors:
			target, score = best, abs(obs - anchors[best])
		for a, anchor_med in anchors.items():
			if a == best:
				continue
			delta, _n = deltas.get((a, best), (None, 0))
			if delta is not None:
				d = abs(obs - (anchor_med - delta))
			elif target is not None:

				continue
			else:
				d = abs(obs - anchor_med) + 1.0
			if score is None or d < score:
				target, score = a, d
		out[qname] = ("attached", target, best)
	return out


def write_fragment_table(path, cluster_id, fragment_results, classified,
						 anchors, deltas, merged_anchors):

	with open(path, "w") as f:
		f.write(
			"Cluster_id\tFragment\tTaxid\tMatched_bases\tIdentity\t"
			"Best_marker\tIs_best\tIs_winner\tStatus\tCredited_to\t"
			"Taxid_is_anchor\tAnchor_level\tDelta_from_credited\t"
			"Predicted_if_own\tPredicted_if_spillover\n"
		)
		for qname in sorted(fragment_results):
			res = fragment_results[qname]
			status, credited, best = classified.get(qname, ("?", None, None))
			winners = set(res["winners"])
			for taxid in sorted(res.get("identity", {})):
				ident = res["identity"][taxid]
				own = anchors.get(taxid)
				delta, _n = deltas.get((credited, taxid), (None, 0)) \
					if credited else (None, 0)
				spill = (anchors[credited] - delta) \
					if (credited in anchors and delta is not None) else None
				f.write(
					f"{cluster_id}\t{qname}\t{taxid}\t"
					f"{res['scores'].get(taxid, '')}\t{round(ident * 100, 3)}\t"
					f"{res['best_marker'].get(taxid, '')}\t"
					f"{'yes' if taxid == best else 'no'}\t"
					f"{'yes' if taxid in winners else 'no'}\t"
					f"{status}\t{credited or 'none'}\t"
					f"{'yes' if taxid in anchors else ('merged:' + merged_anchors[taxid] if taxid in merged_anchors else 'no')}\t"
					f"{round(own * 100, 3) if own is not None else 'NA'}\t"
					f"{round(delta * 100, 3) if delta is not None else 'NA'}\t"
					f"{round(own * 100, 3) if own is not None else 'NA'}\t"
					f"{round(spill * 100, 3) if spill is not None else 'NA'}\n"
				)


def write_pair_table(path, cluster_id, anchors, deltas, merged_anchors,
					 min_separation):
	"""Anchor levels and every measured divergence, with the bridge count behind it."""
	with open(path, "w") as f:
		f.write("Cluster_id\tRecord\tTaxid\tAgainst\tValue\tBridge_fragments\tNote\n")
		for t, lvl in sorted(anchors.items()):
			f.write(f"{cluster_id}\tanchor_level\t{t}\t\t{round(lvl * 100, 3)}\t\t\n")
		for weak, strong in sorted(merged_anchors.items()):
			f.write(f"{cluster_id}\tmerged\t{weak}\t{strong}\t\t\t"
					f"below the {min_separation * 100:.2f} point separation limit\n")
		for (a, c), (d, nb) in sorted(deltas.items()):
			f.write(f"{cluster_id}\tdelta\t{c}\t{a}\t"
					f"{round(d * 100, 3) if d is not None else 'NA'}\t{nb}\t"
					f"{'' if nb >= 3 else 'too few bridge fragments to measure'}\n")


def summarize_evidence(fragment_results, marker_busco,
					   anchor_identity=DEFAULT_ANCHOR_IDENTITY,
					   min_separation=DEFAULT_MIN_SEPARATION,
					   min_anchor_reads=DEFAULT_MIN_ANCHOR_READS,
					   min_bridge_reads=DEFAULT_MIN_BRIDGE_READS,
					   anchor_quantile=DEFAULT_ANCHOR_QUANTILE,
					   anchor_slack=DEFAULT_ANCHOR_SLACK,
					   cluster_id="?", debug_dir=None):

	estimated, level, n_used = estimate_anchor_threshold(
		fragment_results, anchor_quantile, anchor_slack, anchor_identity
	)
	if level is not None:
		floor_binds = estimated <= anchor_identity

		if floor_binds and anchor_identity >= level:
			logger.warning(
				f"{cluster_id}: the anchor_identity floor of "
				f"{anchor_identity:.0%} is at or above the within-species "
				f"level estimated from this cluster's reads ({level:.1%}), so "
				f"the estimate is being discarded and the threshold is "
				f"{estimated:.1%}. Only "
				f"{sum(1 for r in fragment_results.values() if len(r['winners']) == 1 and r.get('identity', {}).get(r['winners'][0], 0) >= estimated)}"
				f" of {n_used} fragments can qualify. anchor_identity is a "
				f"FLOOR now, not the threshold: lower it (0.97 is the default) "
				f"or remove it from the config to let the estimate apply."
			)
		else:
			logger.info(
				f"{cluster_id}: within-species identity estimated at "
				f"{level:.1%} from {n_used} fragments (q{anchor_quantile:g}); "
				f"anchor threshold {estimated:.1%}"
				+ (f", raised to the {anchor_identity:.0%} floor"
				   if floor_binds else "")
			)
	anchor_identity = estimated

	anchors = find_anchors(fragment_results, anchor_identity, min_anchor_reads)

	merged_anchors = {}
	if len(anchors) > 1:
		support = defaultdict(int)
		for res in fragment_results.values():
			if len(res["winners"]) == 1:
				support[res["winners"][0]] += 1
		ranked = sorted(anchors, key=lambda t: (-support[t], t))
		for i, strong in enumerate(ranked):
			if strong in merged_anchors:
				continue
			for weak in ranked[i + 1:]:
				if weak in merged_anchors:
					continue
				d, nb = estimate_delta(
					fragment_results, strong, weak, min_bridge_reads,
					anchor_identity
				)
				if d is not None and abs(d) < min_separation:
					merged_anchors[weak] = strong
					logger.info(
						f"{cluster_id}: {weak} sits {abs(d) * 100:.2f} points "
						f"from {strong}, below the separation limit; their "
						f"reads cannot be told apart."
					)
		for weak in merged_anchors:
			anchors.pop(weak, None)

	candidates = set()
	for res in fragment_results.values():
		candidates.update(res.get("identity", {}))

	deltas = {}
	for a in anchors:
		for c in candidates:
			if c == a:
				continue
			deltas[(a, c)] = estimate_delta(
				fragment_results, a, c, min_bridge_reads, anchor_identity
			)

	classified = classify_fragments(
		fragment_results, anchors, deltas, min_separation, merged_anchors,
		anchor_identity
	)

	ev = defaultdict(lambda: {
		"unique_reads": 0,
		"attached_reads": 0,
		"spillover_reads": 0,
		"undecidable_reads": 0,
		"shared_reads": 0.0,
		"unique_buscos": set(),
		"shared_buscos": set(),
		"identity_sum": 0.0,
		"is_anchor": False,
		"nearest_reference": False,
		"delta_to_anchor": None,
		"bridge_reads": 0,
	})

	undecidable_total = 0
	for qname, (status, credited, best) in classified.items():
		res = fragment_results[qname]
		idents = res.get("identity", {})

		if status == "tied":
			for t in res["winners"]:
				ev[t]["shared_reads"] += 1.0 / len(res["winners"])
				ev[t]["shared_buscos"].add(marker_busco(res["best_marker"][t]))
			continue

		if status == "undecidable":
			undecidable_total += 1
			if best is not None:
				ev[best]["undecidable_reads"] += 1
			continue

		if status == "attached":

			ev[credited]["attached_reads"] += 1
			if best != credited:
				ev[best]["spillover_reads"] += 1
			continue

		# anchor, resolved, unanchored, unverified all support `credited`
		t = credited
		if t not in res["best_marker"]:
			t = best
		ev[t]["unique_reads"] += 1
		ev[t]["identity_sum"] += idents.get(t, 0.0)
		ev[t]["unique_buscos"].add(marker_busco(res["best_marker"][t]))
		if status == "unanchored" and best != t:
			ev[best]["spillover_reads"] += 1
		if status == "unanchored":
			ev[t]["nearest_reference"] = True

	for t in ev:
		ev[t]["unique_buscos"].discard("Collapsed")
		ev[t]["unique_buscos"].discard("Unknown")
		ev[t]["shared_buscos"].discard("Collapsed")
		ev[t]["shared_buscos"].discard("Unknown")
		n = ev[t]["unique_reads"]
		ev[t]["mean_identity"] = (ev[t]["identity_sum"] / n) if n else 0.0
		ev[t]["is_anchor"] = t in anchors
		best_pair = None
		for (a, c), (d, nb) in deltas.items():
			if c == t and d is not None and (best_pair is None or d < best_pair[0]):
				best_pair = (d, nb)
		if best_pair:
			ev[t]["delta_to_anchor"] = round(best_pair[0] * 100, 2)
			ev[t]["bridge_reads"] = best_pair[1]

	by_best = defaultdict(list)
	for res in fragment_results.values():
		if len(res["winners"]) == 1:
			t = res["winners"][0]
			i = res.get("identity", {}).get(t)
			if i is not None:
				by_best[t].append(i)

	def _pct(vals, q):
		if not vals:
			return None
		v = sorted(vals)
		return round(v[min(int(q * len(v)), len(v) - 1)] * 100, 2)

	identity_profile = {
		t: {
			"n": len(v),
			"n_at_anchor_identity": sum(1 for x in v if x >= anchor_identity),
			"p10": _pct(v, 0.10), "p50": _pct(v, 0.50), "p90": _pct(v, 0.90),
		}
		for t, v in by_best.items()
	}

	diagnostics = {
		"anchor_threshold": round(anchor_identity * 100, 3),
		"estimated_level": round(level * 100, 3) if level is not None else None,
		"anchors": {t: round(v * 100, 3) for t, v in anchors.items()},
		"merged_anchors": dict(merged_anchors),
		"undecidable": undecidable_total,
		"identity_profile": identity_profile,
		"deltas": {f"{a}->{c}": (round(d * 100, 2) if d is not None else None, nb)
				   for (a, c), (d, nb) in deltas.items()},
	}

	if debug_dir:
		try:
			os.makedirs(debug_dir, exist_ok=True)
			write_fragment_table(
				os.path.join(debug_dir, f"{cluster_id}_fragments.tsv"),
				cluster_id, fragment_results, classified, anchors, deltas,
				merged_anchors,
			)
			write_pair_table(
				os.path.join(debug_dir, f"{cluster_id}_anchors_and_deltas.tsv"),
				cluster_id, anchors, deltas, merged_anchors, min_separation,
			)
		except Exception as e:
			logger.warning(f"Could not write debug tables for {cluster_id}: {e}")

	dominant_reported = next(
		(t for t, e in ev.items() if e.get("nearest_reference")), None
	)
	diagnostics["nearest_reference"] = dominant_reported

	if anchors:
		attached = sum(e["attached_reads"] for e in ev.values())
		spill = sum(e["spillover_reads"] for e in ev.values())
		logger.info(
			f"{cluster_id}: anchors {sorted(anchors)}; "
			f"{spill} fragments reattached to an anchor, "
			f"{undecidable_total} undecidable"
		)
	else:

		best_seen = max(
			(i for res in fragment_results.values()
			 for i in res.get("identity", {}).values()), default=0.0
		)
		logger.warning(
			f"{cluster_id}: no species reached the anchor identity of "
			f"{anchor_identity:.0%} (best fragment identity seen: "
			f"{best_seen:.1%}). Nothing here can be verified, so the group is "
			f"reported as a nearest-reference match to "
			f"{dominant_reported or 'its best-supported member'} rather than as "
			f"a species identification."
		)

	return dict(ev), diagnostics
