#!/usr/bin/env python
"""
calibrate.py - Decide whether a read stays where it landed or attaches to an anchor.

The question this module exists to answer is the one the previous scoring never
actually asked. When a read's best match is a species other than the abundant
one, is it a read from that species, or is it a read from the abundant species
that had nowhere else to go?

Scoring by "which reference did this read align to best" cannot answer it,
because marker sets are not parallel. After collapsing, an abundant species can
retain a short species-specific fragment of a BUSCO while a relative keeps a
near full-length copy. A read from the abundant species falling outside its own
trimmed marker has only the relative's marker available, wins there by default,
and no comparison ever takes place. That default win is where false positives
come from, and no threshold on alignment score can distinguish it from a real
one, because there is nothing to compare against.

The information needed is already present, in the reads that DO have somewhere
to compare:

  bridge fragment  aligns to markers of two or more species in the group
  orphan fragment  aligns to markers of only one

Bridge fragments measure how far this sample's organism falls when it lands on
each other reference. That measured drop, delta, converts an orphan fragment
from an untestable observation into a testable one: if it is spillover from the
anchor it should sit at (anchor identity - delta); if the candidate is really
present it should sit at the anchor's own within-species identity. Those are two
concrete predictions and the observation falls near one of them.

Two properties make this work where read-count statistics do not:

  * delta is estimated from the ANCHOR's fragments, which are abundant, and
    applied to the CANDIDATE's, which are not. Statistical strength is borrowed
    from the organism that has depth to spare. A handful of candidate reads is
    enough because the two predictions are separated by delta, which is large
    compared to sequencing noise.
  * the within-species prediction comes from the anchor's own median identity in
    this sample, so it absorbs the actual error rate and the actual divergence
    of the sampled strain from its reference. There is no fixed identity cutoff
    to tune.

When delta is too small the two predictions overlap and the fragment is honestly
undecidable, which is reported rather than guessed at.
"""

import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

# A read drawn from the exact reference sequence matches it up to sequencing
# error. This is a property of the instrument and read length, not of taxonomy,
# which is why it is the only fixed constant here.
DEFAULT_ANCHOR_IDENTITY = 0.99

# Minimum separation between the two predictions for a fragment to be decidable.
# The separation IS delta, so this is really "how similar may two species be
# before their reads become indistinguishable". Sequencing noise on a single
# fragment is a few tenths of a point, so one point is comfortably above it.
DEFAULT_MIN_SEPARATION = 0.01

# Matches the presence gate. Two fragments was enough to declare an anchor,
# which then became an attachment target and a calibration yardstick off almost
# no evidence.
DEFAULT_MIN_ANCHOR_READS = 4
DEFAULT_MIN_BRIDGE_READS = 3


def _median(values):
	vals = sorted(values)
	n = len(vals)
	if n == 0:
		return None
	mid = n // 2
	return vals[mid] if n % 2 else (vals[mid - 1] + vals[mid]) / 2.0


def find_anchors(fragment_results, anchor_identity=DEFAULT_ANCHOR_IDENTITY,
				 min_anchor_reads=DEFAULT_MIN_ANCHOR_READS):
	"""
	Species whose own best-matching fragments sit at within-species identity.

	An anchor is an organism that essentially matches its reference. The
	yardstick it provides is what "a read from the species it came from" looks
	like in this sample, with this library, at this strain's divergence.

	Detection counts fragments at within-species identity rather than testing
	the median of all of them. Marker sets contain repetitive and low-complexity
	regions that attract large numbers of poorly-matching reads: in one real
	sample a single Phytophthora infestans marker drew 27415 reads at 60%
	identity, enough to pull that species' median to 93% even though it was
	unambiguously present. A median test throws away obvious anchors whenever
	junk markers outnumber good ones. A count test does not, and the level is
	then taken from the qualifying fragments alone so the junk cannot drag the
	yardstick down either.

	Returns {taxid: median_identity_of_qualifying_fragments}.
	"""
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
	"""
	Measure how far this sample's organism falls when it lands on `candidate`.

	Uses only bridge fragments, those aligning to markers of both species, since
	they are the only ones where both identities are observed on the same
	molecule. Paired that way the comparison is free of any assumption about
	marker coverage.

	Returns (delta, n_bridge). delta is None when there are too few bridge
	fragments to measure it, which is itself a reportable state: the attachment
	cannot be verified from this sample.
	"""
	# Only fragments that are confidently FROM the anchor are usable for
	# calibration. The question being measured is "a read we know came from the
	# anchor -- how much identity does it lose on the candidate's marker", so a
	# fragment that matches the anchor poorly says nothing about that drop and
	# only adds noise. This also keeps junk-marker reads out of the estimate.
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
	"""
	Decide, for every fragment, which taxon it is evidence for.

	The governing rule is that reaching within-species identity is NECESSARY for
	a fragment to be evidence that a non-anchor species is present. Beating the
	spillover prediction is not sufficient on its own, for two reasons found on
	real data:

	  * divergence is not uniform across loci. A single median drop per species
	    pair predicts badly for individual markers, because conserved regions of
	    a BUSCO diverge far less than variable ones. Two Botrytis relatives were
	    called present off reads sitting 2.7 points below the anchor when the
	    pair's median drop was 7 points -- they beat the spillover prediction
	    while being nowhere near within-species identity.
	  * delta cannot always be measured. Four Phytophthora and Tetrahymena
	    relatives had zero bridge fragments, so there was no prediction to beat
	    and their reads were credited by default at 89-91% identity against
	    anchors sitting at 99.8%.

	The measured divergence still does real work: it decides WHICH anchor best
	explains a spillover fragment, and it detects reference pairs too close to
	tell apart. It is simply not the thing that licenses a presence call.

	Statuses:
	  anchor       best match is an anchor; the fragment supports it
	  resolved     best match is not an anchor but reaches within-species
	               identity, so that species really is present too
	  attached     credited to the anchor that best explains it
	  undecidable  the candidate is indistinguishable from a stronger anchor
	  unanchored   no anchor exists in this group, so there is nothing to test
	               against and nothing to attach to

	Returns {qname: (status, credited_taxid, best_taxid)}.
	"""
	merged = merged or {}

	# With no anchor, nothing in this group matches its reference well enough to
	# be called at species level, and there is nothing to attach to. Reads from
	# one absent organism then scatter across whichever references happen to be
	# nearest, and crediting each fragment to its own best match reports one
	# organism as several species -- an organism equidistant from three
	# references produced three calls of 115, 98 and 71 reads. All fragments are
	# instead credited to the single best-supported member, which is reported as
	# a nearest-reference match rather than a species identification.
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

		# The identity requirement applies to anchors too. Being an anchor says
		# the species is present; it says nothing about whether any particular
		# fragment came from it. Exempting a species' fragments once it had
		# cleared the anchor test let a handful of good reads carry hundreds of
		# poor ones: one Blastocystis subtype was reported on 368 fragments
		# averaging 91.5% identity because two of them had reached 99%.
		if obs >= anchor_identity:
			out[qname] = ("anchor" if best in anchors else "resolved", best, best)
			continue

		if not anchors:
			# Credited to the group's best-supported member unless the fragment
			# did not align there at all, in which case it keeps its own best.
			target = dominant if (dominant and dominant in idents) else best
			out[qname] = ("unanchored", target, best)
			continue

		# Below within-species identity, so this fragment is not evidence that
		# `best` is present. Where it is CREDITED is a separate question from
		# what it proves.
		#
		# If `best` is already an anchor, the species is established as present
		# by its own high-identity reads, and this fragment's best explanation
		# is still that species. Giving it to a different anchor would be
		# perverse: in a four-subtype Blastocystis mixture every subtype was an
		# anchor, yet each was handing ~95% of its own fragments to one
		# arbitrary neighbour, which ended up with more reads than the alignment
		# ever gave it. Anchor status governs presence; it must not also decide
		# attribution away from the obvious owner.
		if best in anchors:
			out[qname] = ("attached", best, best)
			continue

		# `best` is not an anchor, so the fragment is spillover. Credit it to
		# whichever anchor explains it best: the one whose measured spillover
		# prediction it sits closest to, falling back to the nearest anchor by
		# identity when no delta could be measured.
		target, score = None, None
		for a, anchor_med in anchors.items():
			delta, _n = deltas.get((a, best), (None, 0))
			if delta is not None:
				d = abs(obs - (anchor_med - delta))
			else:
				d = abs(obs - anchor_med) + 1.0   # unmeasured, deprioritised
			if score is None or d < score:
				target, score = a, d
		out[qname] = ("attached", target, best)
	return out


def summarize_evidence(fragment_results, marker_busco,
					   anchor_identity=DEFAULT_ANCHOR_IDENTITY,
					   min_separation=DEFAULT_MIN_SEPARATION,
					   min_anchor_reads=DEFAULT_MIN_ANCHOR_READS,
					   min_bridge_reads=DEFAULT_MIN_BRIDGE_READS,
					   cluster_id="?"):
	"""
	Full calibration pass over one competition group.

	Replaces scoring by best alignment plus an identity cutoff. Evidence for a
	species is now fragments that were tested against a measured alternative and
	came out on the side of that species being present.

	Returns (evidence, diagnostics).
	"""
	anchors = find_anchors(fragment_results, anchor_identity, min_anchor_reads)

	# The separation test has to apply between anchors as well. Two references
	# only a fraction of a point apart both look like anchors -- reads from
	# either matches both at within-species identity -- but nothing in the data
	# distinguishes them. Keeping both would call two species off one organism,
	# which is the same overcalling in a different disguise. Keep the
	# better-supported anchor and mark the other indistinguishable.
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
			# Credited to the anchor that explains it, and recorded against the
			# reference it landed on so the reassignment stays visible.
			# Spillover counts only fragments actually given AWAY: a fragment
			# that stayed with its own anchor was not reassigned and would
			# otherwise be double-reported as both attached and spilled.
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

	# Per-species identity distribution over the fragments that best-matched it.
	# The spread is what tells a real detection from a species collecting
	# strays, so it is reported rather than reduced to a single mean.
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
		"anchors": {t: round(v * 100, 3) for t, v in anchors.items()},
		"merged_anchors": dict(merged_anchors),
		"undecidable": undecidable_total,
		"identity_profile": identity_profile,
		"deltas": {f"{a}->{c}": (round(d * 100, 2) if d is not None else None, nb)
				   for (a, c), (d, nb) in deltas.items()},
	}

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
		# Loud, because it means the whole verification was inert for this
		# group. If it happens for every group something upstream is wrong --
		# identities not being computed, NM tags missing, or the threshold set
		# above what the data can reach.
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
