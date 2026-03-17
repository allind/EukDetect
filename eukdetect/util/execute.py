from pathlib import Path
from datetime import datetime
from typing import Optional

import logging
import subprocess
import tempfile
import yaml

logger = logging.getLogger(__name__)

class SnakemakeExecutor:
	def __init__(self, config_dict: dict, mode: str = "runall", cores: int = 1, force: bool = False, dry_run: bool = False):
		self.config_dict = config_dict
		self.mode = mode
		self.cores = cores
		self.force = force
		self.dry_run = dry_run
		self.snakefile = self._get_snakefile_path()
	
	def run(self) -> bool:
		#returns true if successful false if fails

		with tempfile.NamedTemporaryFile(
			mode='w',
			suffix='.yml',
			delete=False,
		) as tmp_config:
			yaml.dump(self.config_dict, tmp_config, default_flow_style=False)
			tmp_config_path = tmp_config.name
		
		try:
			snakemake_args = self._build_snakemake_command(tmp_config_path)
			
			# Create logs directory in output folder
			output_dir = Path(self.config_dict["output_dir"])
			log_dir = output_dir / "logs"
			log_dir.mkdir(parents=True, exist_ok=True)
			
			timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
			
			# Include sample name in log filename if single sample
			samples = self.config_dict.get("samples", {})
			if len(samples) == 1:
				sample_name = list(samples.keys())[0]
				log_file = log_dir / f"snakemake_{sample_name}_{timestamp}.log"
			else:
				log_file = log_dir / f"snakemake_{timestamp}.log"
			
			logger.info(f"Running: {' '.join(snakemake_args)}")
			logger.info(f"Snakemake output: {log_file}")
			
			with open(log_file, "w") as log:
				result = subprocess.run(
					snakemake_args,
					stdout=log,
					stderr=subprocess.STDOUT,
					text=True,
				)
			if self.mode == "printaln" and result.returncode == 0:
				self._fix_printaln_output()

			if result.returncode != 0:
				logger.error(f"Snakemake failed with exit code {result.returncode}")
				logger.error(f"Check log file for details: {log_file}")
				return False
			
			return True
		
		finally:
			Path(tmp_config_path).unlink(missing_ok=True)
	
	def _build_snakemake_command(self, config_path: str) -> list:
		
		args = [
			"snakemake",
			"--snakefile", str(self.snakefile),
			"--configfile", config_path,
			"--cores", str(self.cores),
			# Each sample writes to its own output_dir, so pointing Snakemake's
			# working directory there means concurrent jobs each lock a unique
			# directory and never collide with one another.
			"--directory", self.config_dict["output_dir"],
		]
		
		# Add target rule based on mode
		if self.mode == "aln":
			args.append("aln")
		elif self.mode == "analyze":
			args.append("analyze")
			args.extend(["--rerun-triggers", "mtime"])
		elif self.mode == "printaln":
			args.append("printaln")
		# 'all' is default, no need to specify

		# --forceall reruns every rule in the DAG, including runaln.
		# In analyze mode, --force means "overwrite existing filter outputs",
		# which validate.py already handles by deleting them before Snakemake
		# runs. Passing --forceall would unnecessarily re-run alignment.
		# In aln/all modes --forceall is appropriate.
		if self.force and self.mode != "analyze":
			args.append("--forceall")
		
		if self.dry_run:
			args.extend(["--dry-run", "--printshellcmds"])
		
		return args
	
	def _get_snakefile_path(self) -> Path:
		try:
			import eukdetect
			pkg_path = Path(eukdetect.__file__).parent
			snakefile = pkg_path / "rules" / "eukdetect.rules"
			if snakefile.exists():
				return snakefile
		except Exception:
			pass


		snakefile = Path(__file__).parent.parent / "rules" / "eukdetect.rules"
		if snakefile.exists():
			return snakefile


		if "eukdetect_dir" in self.config_dict:
			snakefile = Path(self.config_dict["eukdetect_dir"]) / "rules" / "eukdetect.rules"
			if snakefile.exists():
				return snakefile

		raise FileNotFoundError(
			"Could not locate Snakefile. Make sure EukDetect is properly installed."
		)
	
	def _fix_printaln_output(self) -> None:
		output_dir = Path(self.config_dict["output_dir"])
		cmd_file = output_dir / "alignment_commands.txt"
		
		if not cmd_file.exists():
			return
		
		logger.debug("Fixing printaln output formatting")
		
		with open(cmd_file) as f:
			lines = f.readlines()
		
		fixed_lines = [line.replace("{}", "'") for line in lines]
		
		with open(cmd_file, "w") as f:
			f.writelines(fixed_lines)
