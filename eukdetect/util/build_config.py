
from Bio import SeqIO
from pathlib import Path
from typing import Dict, Optional

import logging
import gzip

#set up logging
logger = logging.getLogger(__name__)

class ConfigBuilder:

	def __init__(self, samples: Dict[str,dict], output_dir: str, database_dir: str, database_prefix: str = "eukdb_clust97", paired_end: bool = True, readlen: Optional[int] = None, bowtie2_cores: int = 1):
		self.samples = samples
		self.output_dir = output_dir
		self.database_dir = database_dir
		self.database_prefix = database_prefix
		self.paired_end = paired_end
		self.readlen = readlen
		self.bowtie2_cores = bowtie2_cores

	def build(self) -> dict:

		if self.readlen is None:
			self.readlen = self._detect_readlen()
			logger.info(f"Auto-detected read length: {self.readlen} bp")

		fwd_suffix, rev_suffix, se_suffix = self._determine_suffixes()

		fastq_dir = self._get_fastq_dir()
		eukdetect_dir = self._get_eukdetect_dir()

		config = {
			"output_dir": str(Path(self.output_dir).absolute()),
			"paired_end": self.paired_end,
			"fwd_suffix": fwd_suffix,
			"rev_suffix": rev_suffix,
			"se_suffix": se_suffix,
			"readlen": self.readlen,
			"fq_dir": fastq_dir,
			"database_dir": str(Path(self.database_dir).absolute()),
			"database_prefix": self.database_prefix,
			"eukdetect_dir": eukdetect_dir,
			"bowtie2_cores": self.bowtie2_cores,
			"samples": {name: None for name in self.samples.keys()}
		}

		return config

	def _detect_readlen(self) -> int:
		first_sample = list(self.samples.values())[0]
		fastq_path = first_sample["reads1"]
		logger.debug(f"Detecting read length from {fastq_path}")
		
		counter = 0
		bases = 0
		max_reads = 10000
		
		try:
			if fastq_path.endswith(".gz"):
				with gzip.open(fastq_path, "rt") as handle:
					for record in SeqIO.parse(handle, "fastq"):
						if counter >= max_reads:
							break
						counter += 1
						bases += len(record.seq)
			else:
				for record in SeqIO.parse(fastq_path, "fastq"):
					if counter >= max_reads:
						break
					counter += 1
					bases += len(record.seq)
			
			if counter == 0:
				raise ValueError("Empty FASTQ")

			return int(bases / counter)

		except ValueError:
			logger.error(f"Input FASTQ appears to be empty: {fastq_path}")
			logger.error("Specify read length manually with --readlen to bypass auto-detection.")
			raise
		except (OSError, IOError) as e:
			logger.warning(f"Could not auto-detect read length: {e}")
			logger.warning("Defaulting to 150 bp. Use --readlen to specify manually.")
			return 150


	def _determine_suffixes(self) -> tuple:

		
		first_sample = list(self.samples.values())[0]
		first_name = list(self.samples.keys())[0]
		fastq1 = Path(first_sample["reads1"]).name
		
		if self.paired_end and "reads2" in first_sample:
			fastq2 = Path(first_sample["reads2"]).name
			
			# Remove sample name to get suffix
			fwd_suffix = fastq1.replace(first_name, "", 1)
			rev_suffix = fastq2.replace(first_name, "", 1)
			se_suffix = ".fastq.gz"  # default
			
			logger.debug(f"Detected paired-end suffixes: {fwd_suffix}, {rev_suffix}")
		else:
			# Single-end
			se_suffix = fastq1.replace(first_name, "", 1)
			fwd_suffix = "_1.fastq.gz"  # defaults
			rev_suffix = "_2.fastq.gz"
			
			logger.debug(f"Detected single-end suffix: {se_suffix}")
		
		return fwd_suffix, rev_suffix, se_suffix


	def _get_fastq_dir(self) -> str:

		first_sample = list(self.samples.values())[0]
		fastq_path = Path(first_sample["reads1"])
		return str(fastq_path.parent.absolute())
	
	def _get_eukdetect_dir(self) -> str:
		
		try:
			import eukdetect
			pkg_path = Path(eukdetect.__file__).parent
			return str(pkg_path.absolute())
		except Exception:
			# Fallback for development
			logger.warning("Could not determine package installation path")
			return str(Path(__file__).parent.parent.absolute())


