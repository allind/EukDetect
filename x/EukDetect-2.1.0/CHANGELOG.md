# Changelog

## 2.1.0

### Output format

These are breaking changes to file names and column contents. Scripts written
against 2.0 output need updating.

- Tab-separated output files now end in `.tsv` instead of `.txt`. Affected:
  `{sample}_filtered_hits_table`, `{sample}_filtered_hits_eukfrac`,
  `filtering/{sample}_all_hits_table`, `filtering/{sample}_read_identities`,
  `filtering/{sample}_reassignment_report`,
  `filtering/{sample}_read_counts_and_mismatches`. `alignment_commands.txt` is
  not tabular and keeps its `.txt` extension. Database file names are unchanged,
  so existing databases work without modification.
- Missing values are written as `NA` in every column of every file. Previously
  `Reassigned_genomes` used `None` and the `AS` column of the read identity
  table used an empty string. `None` and `NA` are both null tokens to
  `pandas.read_csv`, so a `Reassigned_genomes` column holding both was silently
  read as all-missing and the genome accessions were lost.
- Percentage columns are written as bare numbers on a 0-100 scale. `PID_aligned`
  and `PID_reassigned` in the filtered hits table, and `Total_marker_coverage`
  and `Percent_identity` in the all-hits table, previously carried a `%` suffix
  that forced them to parse as strings.
- `PID_reassigned` is `NA` when no reads were reassigned to a taxon. It was
  previously `0.0%`, which averaged into summaries as a real measurement.

### Packaging

- Version is defined once, in `eukdetect/__init__.py`, and read from there by
  the build. `setup.py`, `pyproject.toml` and `meta.yaml` previously declared
  three different versions.
- `setup.py` removed. `pyproject.toml` is the only build configuration.
- Runtime dependencies are declared in `pyproject.toml`, so `pip install
  eukdetect` now installs a working environment. `pandas` and `numpy` were
  listed as conda dependencies but are not imported anywhere, and have been
  dropped.
- `eukdetect/envs/eukdetect2_environment.yml` is now included in the built
  package. It was declared as package data for a directory that was not a
  package, so it shipped in neither the wheel nor the sdist.
- `six` added as an explicit dependency. `ete3` imports it without declaring it.
- The root `environment.yaml` has been removed. It was an export of a
  development environment named `eukdetect-realign`, pinned to exact build
  strings, and referenced nowhere. The conda environment used by the README's
  source install, `eukdetect/envs/eukdetect2_environment.yml`, is unchanged.

### Python 3.13

The `ete3` package on PyPI imports the standard library `cgi` module, which was
removed in Python 3.13. A `pip install` on 3.13 therefore builds and installs
cleanly and then fails at the first taxonomy step. The conda-forge build of
`ete3` does not have this problem, so conda installs on 3.13 work.

Rather than cap `requires-python` and break conda installs on 3.13, EukDetect
now catches this specific failure and exits with an instruction to use conda or
Python 3.12. Classifiers list support through 3.12.

### Documentation

- Every output file and column is documented in the README with its type and
  its missing-value behaviour.
- The normalized table's relative abundance column is documented as
  `{sample}_EukFrac`, the name the code has always written. The README
  previously called it `{sample}_RelEuk`.
- The column list for `{sample}_filtered_hits_table` described columns belonging
  to the all-hits table.
