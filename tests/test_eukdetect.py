"""
Unit tests for EukDetect
Run with: pytest tests/test_eukdetect.py --basetemp=./tmp_test -v
"""

import pytest
from pathlib import Path
import sys
import argparse
import yaml
import tempfile
import os

# Add parent directory to path to import eukdetect modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from eukdetect.wrappers.runall import validate_cores, validate_sample_name, _parse_samples
from eukdetect.util.build_config import ConfigBuilder
from eukdetect.util.validate import check_fastq_files, check_database, validate_inputs
from eukdetect.util.execute import SnakemakeExecutor


# Test fixtures - create dummy FASTQ files
@pytest.fixture
def test_fastq_files(tmp_path):
	"""Create dummy FASTQ files for testing"""
	# Single-end file
	se_file = tmp_path / "test.fastq.gz"
	se_file.write_bytes(b'\x1f\x8b\x08\x00')  # Minimal gzip header
	
	# Paired-end files
	pe_file1 = tmp_path / "test_1.fastq.gz"
	pe_file2 = tmp_path / "test_2.fastq.gz"
	pe_file1.write_bytes(b'\x1f\x8b\x08\x00')
	pe_file2.write_bytes(b'\x1f\x8b\x08\x00')
	
	return {
		'se': str(se_file),
		'pe1': str(pe_file1),
		'pe2': str(pe_file2),
		'dir': str(tmp_path)
	}


@pytest.fixture
def test_database(tmp_path):
	"""Create minimal dummy database files"""
	db_dir = tmp_path / "database"
	db_dir.mkdir()
	
	# Create required database files
	prefix = "eukdb"
	
	# Bowtie2 index files (.bt2l variant)
	for suffix in ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']:
		(db_dir / f"{prefix}{suffix}").touch()
	
	# Taxa database
	(db_dir / "taxa.sqlite").touch()
	
	# Taxon list
	(db_dir / "taxon_list.txt").write_text("taxon1\ntaxon2\n")
	
	# Length file
	(db_dir / f"{prefix}_lengths.txt").write_text("seq1\t1000\nseq2\t2000\n")
	
	# Additional required database files
	(db_dir / l").touch()
	(db_dir / "specific_and_inherited_markers_per_taxid.txt").touch()
	(db_dir / "busco_taxid_genome_link.txt").touch()
	(db_dir / "taxid_and_genome_cumulativelength.txt").touch()
	
	return {
		'dir': str(db_dir),
		'prefix': prefix
	}


class TestValidateCores:
	"""Test cores validation function"""
	
	def test_valid_cores(self):
		"""Test that valid core counts are accepted"""
		assert validate_cores("1") == 1
		assert validate_cores("4") == 4
		assert validate_cores("16") == 16
	
	def test_negative_cores_rejected(self):
		"""Test that negative values are rejected"""
		with pytest.raises(argparse.ArgumentTypeError, match="must be at least 1"):
			validate_cores("0")
		
		with pytest.raises(argparse.ArgumentTypeError, match="must be at least 1"):
			validate_cores("-1")
	
	def test_non_integer_rejected(self):
		"""Test that non-integers are rejected"""
		with pytest.raises(argparse.ArgumentTypeError, match="must be an integer"):
			validate_cores("abc")
		
		with pytest.raises(argparse.ArgumentTypeError, match="must be an integer"):
			validate_cores("3.14")
	
	def test_warning_for_excessive_cores(self, caplog):
		"""Test that warning is issued for excessive core counts"""
		import os
		max_cores = os.cpu_count() or 1
		excessive = (max_cores * 2) + 1
		
		result = validate_cores(str(excessive))
		assert result == excessive
		if max_cores:
			assert "performance issues" in caplog.text.lower()


class TestValidateSampleName:
	"""Test sample name validation function"""
	
	def test_valid_sample_names(self):
		"""Test that valid sample names are accepted"""
		assert validate_sample_name("sample1") == "sample1"
		assert validate_sample_name("my_sample_01") == "my_sample_01"
		assert validate_sample_name("sample-2024.01") == "sample-2024.01"
		assert validate_sample_name("ERR4097171") == "ERR4097171"
		assert validate_sample_name("test") == "test"
	
	def test_alphanumeric_only(self):
		"""Test that only alphanumeric characters are allowed"""
		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample with spaces")
		
		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample/path")
		
		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample@email")
		
		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample#1")
	
	def test_no_leading_period_or_dash(self):
		"""Test that names cannot start with . or -"""
		with pytest.raises(ValueError, match="cannot start with"):
			validate_sample_name(".hidden")
		
		with pytest.raises(ValueError, match="cannot start with"):
			validate_sample_name("-sample")
	
	def test_no_path_traversal(self):
		"""Test that path traversal is blocked"""
		# Names with / are caught by the valid characters check first
		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("../../../etc/passwd")

		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample/../other")

		with pytest.raises(ValueError, match="Invalid sample name"):
			validate_sample_name("sample\\windows\\path")


class TestParseSamples:
	"""Test sample parsing from command line and files"""
	
	def test_parse_single_sample_with_name(self):
		"""Test parsing single sample with explicit name"""
		samples = _parse_samples(
			reads1=["data/sample_R1.fq.gz"],
			reads2=["data/sample_R2.fq.gz"],
			sample_name=["my_sample"],
			samples_file=None
		)
		
		assert "my_sample" in samples
		assert samples["my_sample"]["reads1"] == "data/sample_R1.fq.gz"
		assert samples["my_sample"]["reads2"] == "data/sample_R2.fq.gz"
	
	def test_parse_single_sample_auto_name(self):
		"""Test auto-detection of sample name from filename"""
		samples = _parse_samples(
			reads1=["path/to/ERR4097171_R1.fastq.gz"],
			reads2=["path/to/ERR4097171_R2.fastq.gz"],
			sample_name=[],
			samples_file=None
		)
		
		assert "ERR4097171" in samples
		assert samples["ERR4097171"]["reads1"] == "path/to/ERR4097171_R1.fastq.gz"
	
	def test_parse_auto_name_removes_suffixes(self):
		"""Test that common suffixes are removed from auto-detected names"""
		test_cases = [
			("sample_R1.fastq.gz", "sample"),
			("sample_1.fq.gz", "sample"),
			("sample_fwd.fastq", "sample"),
			("sample.fastq.gz", "sample"),
			("test_1.fastq.gz", "test"),
		]
		
		for filename, expected_name in test_cases:
			samples = _parse_samples(
				reads1=[filename],
				reads2=[],
				sample_name=[],
				samples_file=None
			)
			assert expected_name in samples, f"Failed for {filename}"
	
	def test_parse_single_end_sample(self):
		"""Test parsing single-end sample (no R2)"""
		samples = _parse_samples(
			reads1=["test.fastq.gz"],
			reads2=[],
			sample_name=["test"],
			samples_file=None
		)
		
		assert "test" in samples
		assert "reads1" in samples["test"]
		assert "reads2" not in samples["test"]
	
	def test_invalid_sample_name_rejected(self):
		"""Test that invalid sample names are rejected"""
		with pytest.raises(ValueError, match="Invalid sample name"):
			_parse_samples(
				reads1=["data/sample_R1.fq.gz"],
				reads2=["data/sample_R2.fq.gz"],
				sample_name=["invalid/name"],
				samples_file=None
			)
	
	def test_parse_from_file(self, tmp_path):
		"""Test parsing samples from TSV file"""
		samples_file = tmp_path / "samples.tsv"
		samples_file.write_text(
			"sample_name\treads1\treads2\n"
			"sample1\t/path/to/sample1_R1.fq.gz\t/path/to/sample1_R2.fq.gz\n"
			"sample2\t/path/to/sample2_R1.fq.gz\t/path/to/sample2_R2.fq.gz\n"
		)
		
		samples = _parse_samples(
			reads1=[],
			reads2=[],
			sample_name=[],
			samples_file=str(samples_file)
		)
		
		assert len(samples) == 2
		assert "sample1" in samples
		assert "sample2" in samples
		assert samples["sample1"]["reads1"] == "/path/to/sample1_R1.fq.gz"
		assert samples["sample1"]["reads2"] == "/path/to/sample1_R2.fq.gz"
	
	def test_parse_from_file_single_end(self, tmp_path):
		"""Test parsing single-end samples from TSV file"""
		samples_file = tmp_path / "samples.tsv"
		samples_file.write_text(
			"sample_name\treads1\n"
			"test\t/path/to/test.fastq.gz\n"
		)
		
		samples = _parse_samples(
			reads1=[],
			reads2=[],
			sample_name=[],
			samples_file=str(samples_file)
		)
		
		assert "test" in samples
		assert "reads1" in samples["test"]
		assert "reads2" not in samples["test"]
	
	def test_invalid_sample_name_in_file_rejected(self, tmp_path):
		"""Test that invalid sample names in TSV are rejected"""
		samples_file = tmp_path / "samples.tsv"
		samples_file.write_text(
			"sample_name\treads1\treads2\n"
			"../bad/name\t/path/to/R1.fq.gz\t/path/to/R2.fq.gz\n"
		)
		
		with pytest.raises(ValueError, match="Invalid sample name"):
			_parse_samples(
				reads1=[],
				reads2=[],
				sample_name=[],
				samples_file=str(samples_file)
			)


class TestConfigBuilding:
	"""Test configuration building"""
	
	def test_build_config_paired_end(self, test_fastq_files, test_database):
		"""Test building config for paired-end sample"""
		samples = {
			"test": {
				"reads1": test_fastq_files['pe1'],
				"reads2": test_fastq_files['pe2']
			}
		}
		
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		assert config['paired_end'] is True
		assert config['readlen'] == 100
		assert config['bowtie2_cores'] == 8
		assert config['database_prefix'] == test_database['prefix']
		assert 'test' in config['samples']
		assert config['fwd_suffix'] == '_1.fastq.gz'
		assert config['rev_suffix'] == '_2.fastq.gz'
	
	def test_build_config_single_end(self, test_fastq_files, test_database):
		"""Test building config for single-end sample"""
		samples = {
			"test": {
				"reads1": test_fastq_files['se']
			}
		}
		
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=False,
			readlen=100,
			bowtie2_cores=4
		)
		
		config = builder.build()
		
		assert config['paired_end'] is False
		assert config['readlen'] == 100
		assert config['bowtie2_cores'] == 4
		assert 'test' in config['samples']
		assert config['se_suffix'] == '.fastq.gz'
	
	def test_config_has_required_keys(self, test_fastq_files, test_database):
		"""Test that config contains all required keys"""
		samples = {"test": {"reads1": test_fastq_files['pe1'], "reads2": test_fastq_files['pe2']}}
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		required_keys = [
			'output_dir', 'paired_end', 'fwd_suffix', 'rev_suffix', 'se_suffix',
			'readlen', 'fq_dir', 'database_dir', 'database_prefix', 'bowtie2_cores', 'samples'
		]
		
		for key in required_keys:
			assert key in config, f"Missing required key: {key}"


class TestFileValidation:
	"""Test file validation functions"""
	
	def test_check_database_valid(self, test_database):
		"""Test that valid database passes validation"""
		config = {
			'database_dir': test_database['dir'],
			'database_prefix': test_database['prefix']
		}
		
		# Should not raise any exception
		check_database(config)
	
	def test_check_database_missing_files(self, tmp_path):
		"""Test that missing database files are detected"""
		db_dir = tmp_path / "incomplete_db"
		db_dir.mkdir()
		
		# Only create some files, not all
		(db_dir / "eukdb.1.bt2l").touch()
		
		config = {
			'database_dir': str(db_dir),
			'database_prefix': 'eukdb'
		}
		
		with pytest.raises(ValueError, match="Missing database files"):
			check_database(config)
	
	def test_check_fastq_files_paired_end(self, test_fastq_files):
		"""Test FASTQ file validation for paired-end"""
		config = {
			'samples': {'test': None},
			'fq_dir': test_fastq_files['dir'],
			'paired_end': True,
			'fwd_suffix': '_1.fastq.gz',
			'rev_suffix': '_2.fastq.gz'
		}
		
		# Should not raise exception
		check_fastq_files(config)
	
	def test_check_fastq_files_single_end(self, test_fastq_files):
		"""Test FASTQ file validation for single-end"""
		config = {
			'samples': {'test': None},
			'fq_dir': test_fastq_files['dir'],
			'paired_end': False,
			'se_suffix': '.fastq.gz'
		}
		
		# Should not raise exception
		check_fastq_files(config)
	
	def test_check_fastq_files_missing(self, tmp_path):
		"""Test that missing FASTQ files are detected"""
		config = {
			'samples': {'nonexistent': None},
			'fq_dir': str(tmp_path),
			'paired_end': True,
			'fwd_suffix': '_1.fastq.gz',
			'rev_suffix': '_2.fastq.gz'
		}
		
		with pytest.raises(ValueError, match="Missing input fastq files"):
			check_fastq_files(config)
	
	def test_validate_inputs_complete(self, test_fastq_files, test_database):
		"""Test complete input validation"""
		samples = {"test": {"reads1": test_fastq_files['pe1'], "reads2": test_fastq_files['pe2']}}
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		# Should not raise exception
		validate_inputs(config, mode='all', force=False)


class TestSnakemakeExecutor:
	"""Test Snakemake execution preparation"""
	
	def test_executor_initialization(self, test_fastq_files, test_database):
		"""Test that SnakemakeExecutor initializes correctly"""
		samples = {"test": {"reads1": test_fastq_files['pe1'], "reads2": test_fastq_files['pe2']}}
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		executor = SnakemakeExecutor(
			config_dict=config,
			mode='all',
			cores=8,
			force=False,
			dry_run=True
		)
		
		assert executor.config_dict == config
		assert executor.mode == 'all'
		assert executor.cores == 8
		assert executor.force is False
		assert executor.dry_run is True
	
	def test_build_snakemake_command_basic(self, test_fastq_files, test_database):
		"""Test building basic Snakemake command"""
		samples = {"test": {"reads1": test_fastq_files['pe1'], "reads2": test_fastq_files['pe2']}}
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		executor = SnakemakeExecutor(
			config_dict=config,
			mode='all',
			cores=8,
			force=False,
			dry_run=False
		)
		
		# Create a temporary config file
		with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
			yaml.dump(config, f)
			temp_config = f.name
		
		try:
			cmd = executor._build_snakemake_command(temp_config)
			
			assert 'snakemake' in cmd
			assert '--cores' in cmd
			assert '8' in cmd
			assert '--configfile' in cmd
			assert temp_config in cmd
		finally:
			os.unlink(temp_config)
	
	def test_build_snakemake_command_with_mode(self, test_fastq_files, test_database):
		"""Test building Snakemake command with different modes"""
		samples = {"test": {"reads1": test_fastq_files['pe1'], "reads2": test_fastq_files['pe2']}}
		output_dir = Path(test_fastq_files['dir']) / "output"
		output_dir.mkdir()
		
		builder = ConfigBuilder(
			samples=samples,
			output_dir=str(output_dir),
			database_dir=test_database['dir'],
			database_prefix=test_database['prefix'],
			paired_end=True,
			readlen=100,
			bowtie2_cores=8
		)
		
		config = builder.build()
		
		# Test each mode
		for mode in ['all', 'aln', 'analyze', 'printaln']:
			executor = SnakemakeExecutor(
				config_dict=config,
				mode=mode,
				cores=8,
				force=False,
				dry_run=False
			)
			
			with tempfile.NamedTemporaryFile(mode='w', suffix='.yml', delete=False) as f:
				yaml.dump(config, f)
				temp_config = f.name
			
			try:
				cmd = executor._build_snakemake_command(temp_config)
				
				# 'all' mode doesn't add rule name, others do
				if mode != 'all':
					assert mode in cmd
			finally:
				os.unlink(temp_config)


if __name__ == "__main__":
	# Allow running with: python test_eukdetect.py
	pytest.main([__file__, "-v"])
