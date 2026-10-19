#!/usr/bin/env python
"""
make_identity_cache.py - Build <database_dir>/identity_cache.json ahead of time.

This is optional. EukDetect builds the cache itself when it is missing or stale,
so a normal run needs nothing here. Use this to pay the one-off cost up front,
or to write the cache as a user who can write to the database directory when the
users who run analyses cannot.

Usage:
    python make_identity_cache.py --database_dir /path/to/eukdb
    python make_identity_cache.py --database_dir /path/to/eukdb --force
"""

import argparse
import importlib.util
import logging
import os
import sys

logging.basicConfig(
	format="%(asctime)s [%(levelname)s] %(message)s",
	datefmt="%Y-%m-%d %H:%M:%S",
	level=logging.INFO,
)
logger = logging.getLogger(__name__)

TAXID_LINK = "busco_taxid_genome_link.txt"
TAXDB = "taxa.sqlite"


def _load_ani_clusters():
	"""Import eukdetect.util.ani_clusters whether or not the package is installed."""
	try:
		from eukdetect.util import ani_clusters
		return ani_clusters
	except ImportError:
		path = os.path.join(
			os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
			"eukdetect", "util", "ani_clusters.py",
		)
		if not os.path.exists(path):
			logger.error(f"Cannot find ani_clusters.py (looked in {path})")
			sys.exit(1)
		spec = importlib.util.spec_from_file_location("ani_clusters", path)
		mod = importlib.util.module_from_spec(spec)
		spec.loader.exec_module(mod)
		return mod


def main(argv=None):
	p = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
	)
	p.add_argument("--database_dir", required=True)
	p.add_argument("--force", action="store_true",
				   help="Rebuild even if the cache is present and current")
	args = p.parse_args(argv)

	ac = _load_ani_clusters()
	db = args.database_dir
	taxid_link = os.path.join(db, TAXID_LINK)
	taxdb = os.path.join(db, TAXDB)

	for path in (taxid_link, taxdb):
		if not os.path.exists(path):
			logger.error(f"Missing required database file: {path}")
			sys.exit(1)

	if not args.force:
		existing = ac.load_identity_cache(db, taxid_link)
		if existing:
			logger.info(
				f"{os.path.join(db, ac.IDENTITY_CACHE_FILENAME)} is already "
				f"current ({existing.get('n_taxa')} taxa). Use --force to rebuild."
			)
			return

	if not os.access(db, os.W_OK):
		logger.error(f"Database directory is not writable: {db}")
		sys.exit(1)

	from ete3 import NCBITaxa
	logger.info(f"Loading taxonomy from {taxdb}")
	ncbi = NCBITaxa(taxdb)

	payload = ac.build_identity_cache(db, taxid_link, ncbi)
	for rank in payload["cached_ranks"]:
		n = sum(1 for e in payload["ranks"].values() if rank in e)
		logger.info(f"  {rank}: {n} taxa")


if __name__ == "__main__":
	main()
