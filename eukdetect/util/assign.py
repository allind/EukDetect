#!/usr/bin/env python


import logging
from collections import defaultdict

logger = logging.getLogger(__name__)

DEFAULT_MIN_UNIQUE_READS = 4
DEFAULT_MIN_UNIQUE_MARKERS = 2
# Kept in step with realign.DEFAULT_MIN_UNIQUE_IDENTITY. The safe value is the
# default so a caller that forgets to pass it gets filtering rather than none.
DEFAULT_MIN_UNIQUE_IDENTITY = 0.97
EM_MAX_ITERS = 500
EM_TOLERANCE = 1e-8


def apply_presence_gate(evidence, min_unique_reads=DEFAULT_MIN_UNIQUE_READS,
						min_unique_markers=DEFAULT_MIN_UNIQUE_MARKERS):

	present, failed = [], []
	for taxid, ev in evidence.items():

		if (ev["unique_reads"] >= min_unique_reads
				and len(ev["unique_buscos"]) >= min_unique_markers):
			present.append(taxid)
		else:
			failed.append(taxid)
	return sorted(present), sorted(failed)


def run_em(fragment_results, present, marker_lengths,
		   min_identity=DEFAULT_MIN_UNIQUE_IDENTITY,
		   max_iters=EM_MAX_ITERS, tolerance=EM_TOLERANCE):

	present_set = set(present)
	if not present_set:
		return {}, len(fragment_results)

	# Unique fragments anchor the abundances; only ties need iterating.
	anchors = defaultdict(float)
	ties = []
	unassignable = 0

	for res in fragment_results.values():
		idents = res.get("identity", {})
		cands = [
			t for t in res["winners"]
			if t in present_set and idents.get(t, 1.0) >= min_identity
		]
		if not cands:
			unassignable += 1
			continue
		if len(cands) == 1:
			anchors[cands[0]] += 1.0
		else:
			ties.append(cands)

	def eff_len(t):
		return max(marker_lengths.get(t, 0), 1)

	theta = {t: anchors.get(t, 0.0) + 1e-9 for t in present_set}
	total = sum(theta.values())
	theta = {t: v / total for t, v in theta.items()}

	for it in range(max_iters):
		counts = {t: anchors.get(t, 0.0) for t in present_set}
		for cands in ties:
			weights = {t: theta[t] / eff_len(t) for t in cands}
			denom = sum(weights.values())
			if denom <= 0:
				share = 1.0 / len(cands)
				for t in cands:
					counts[t] += share
				continue
			for t in cands:
				counts[t] += weights[t] / denom

		total = sum(counts.values())
		if total <= 0:
			break
		new_theta = {t: counts[t] / total for t in present_set}
		delta = max(abs(new_theta[t] - theta[t]) for t in present_set)
		theta = new_theta
		if delta < tolerance:
			logger.debug(f"EM converged after {it + 1} iterations (delta={delta:.2e})")
			break
	else:
		logger.warning(
			f"EM hit the {max_iters}-iteration cap without converging; "
			f"using the last estimate."
		)

	# Recover final read counts from the converged responsibilities.
	assigned = {t: anchors.get(t, 0.0) for t in present_set}
	for cands in ties:
		weights = {t: theta[t] / eff_len(t) for t in cands}
		denom = sum(weights.values())
		for t in cands:
			assigned[t] += (weights[t] / denom) if denom > 0 else (1.0 / len(cands))

	return assigned, unassignable


def lca_of(taxids, ncbi):

	lineages = []
	for t in taxids:
		try:
			lineages.append(ncbi.get_lineage(int(t)))
		except Exception as e:
			logger.warning(f"No lineage for taxid {t}: {e}")
	if not lineages:
		return None

	shared = lineages[0]
	for lin in lineages[1:]:
		limit = min(len(shared), len(lin))
		common = []
		for i in range(limit):
			if shared[i] == lin[i]:
				common.append(shared[i])
			else:
				break
		shared = common
		if not shared:
			return None
	return str(shared[-1]) if shared else None


def resolve_cluster(cluster_id, member_taxids, fragment_results, evidence,
					marker_lengths, ncbi=None,
					min_unique_reads=DEFAULT_MIN_UNIQUE_READS,
					min_unique_markers=DEFAULT_MIN_UNIQUE_MARKERS,
					min_identity=DEFAULT_MIN_UNIQUE_IDENTITY,
					min_cluster_reads_for_lca=None):

	if not evidence:
		return []

	present, failed = apply_presence_gate(
		evidence, min_unique_reads, min_unique_markers
	)

	if present:
		assigned, unassignable = run_em(
			fragment_results, present, marker_lengths, min_identity=min_identity
		)
		# Fragments the calibration reattached belong to the anchor that
		# explains them, not to the reference they happened to align to.
		for t in present:
			assigned[t] = assigned.get(t, 0.0) + evidence[t].get("attached_reads", 0)
		if unassignable:
			logger.debug(
				f"{cluster_id}: {unassignable} fragments belonged only to species "
				f"that failed the gate; not forced onto a present species."
			)
		calls = []
		for t in present:
			ev = evidence[t]
			calls.append({
				"taxid": t,
				"cluster_id": cluster_id,
				# nearest_reference means nothing in the group matched its
				# reference well enough to be identified to species; this is the
				# closest thing the database has, at the identity reported.
				"resolution": ("nearest_reference"
							   if ev.get("nearest_reference") else "species"),
				"unique_reads": ev["unique_reads"],
				"unique_markers": len(ev["unique_buscos"]),
				"shared_reads": round(ev["shared_reads"], 3),
				"attached_reads": ev.get("attached_reads", 0),
				"spillover_reads": ev.get("spillover_reads", 0),
				"undecidable_reads": ev.get("undecidable_reads", 0),
				"is_anchor": ev.get("is_anchor", False),
				"delta_to_anchor": ev.get("delta_to_anchor"),
				"bridge_reads": ev.get("bridge_reads", 0),
				"mean_identity": round(ev.get("mean_identity", 0.0) * 100, 2),
				"assigned_reads": round(assigned.get(t, 0.0), 3),
				"members": sorted(member_taxids),
			})
		if failed:
			logger.info(
				f"{cluster_id}: {len(present)} species called, "
				f"{len(failed)} lacked unique evidence ({','.join(failed)})"
			)
		return calls

	total_reads = len(fragment_results)
	undecidable = sum(e.get("undecidable_reads", 0) for e in evidence.values())
	if undecidable:
		logger.info(
			f"{cluster_id}: {undecidable} fragments could not distinguish their "
			f"candidates because the measured divergence between them was below "
			f"the separation limit."
		)
	threshold = (min_cluster_reads_for_lca
				 if min_cluster_reads_for_lca is not None
				 else min_unique_reads)
	if total_reads < threshold:
		logger.debug(f"{cluster_id}: {total_reads} fragments, below reporting floor.")
		return []

	if len(member_taxids) == 1 or ncbi is None:
		return []

	observed = sorted(evidence.keys())
	lca = lca_of(observed, ncbi)
	if lca is None:
		logger.warning(f"{cluster_id}: could not compute an LCA; group dropped.")
		return []

	total_shared = sum(ev["shared_reads"] for ev in evidence.values())

	member_lengths = [marker_lengths.get(t, 0) for t in observed]
	member_lengths = [L for L in member_lengths if L > 0]
	mean_length = (sum(member_lengths) / len(member_lengths)) if member_lengths else 0

	logger.info(
		f"{cluster_id}: no member has unique evidence over {total_reads} "
		f"fragments; reporting at LCA taxid {lca} as unresolved_cluster. "
		f"Candidates were {','.join(observed)}."
	)
	return [{
		"taxid": lca,
		"cluster_id": cluster_id,
		"resolution": "unresolved_cluster",
		"unique_reads": 0,
		"unique_markers": 0,
		"shared_reads": round(total_shared, 3),
		"attached_reads": 0,
		"spillover_reads": 0,
		"undecidable_reads": sum(e.get("undecidable_reads", 0) for e in evidence.values()),
		"is_anchor": False,
		"delta_to_anchor": None,
		"bridge_reads": 0,
		"mean_identity": 0.0,
		"assigned_reads": float(total_reads),
		"marker_length": mean_length,
		"members": observed,
	}]
