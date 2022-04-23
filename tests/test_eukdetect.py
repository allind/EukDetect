#!/usr/bin/env python

import unittest
import os
import shutil
import subprocess
import sys
import shlex

def run(command):
	process = subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	return(process.returncode, process.stdout, process.stderr)

class _01_HelpText(unittest.TestCase):

	def test_class(self):
		command = "eukdetect -h"
		code, stdout, stderr = run(command)
		self.assertTrue(code == 0, msg=stdout)


class _02_RunRunall(unittest.TestCase):
	def test_class(self):
		command = 'eukdetect --mode runall --configfile tests/configfile_for_tests.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)

class _03_RunAlnCmd(unittest.TestCase):

	def test_class(self):
		command = 'eukdetect --mode printaln --configfile tests/configfile_for_tests.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)


class _04_RunAln_and_filter(unittest.TestCase):

	def test_class(self):
		command = 'eukdetect --mode aln --configfile tests/configfile_for_tests.yml --force'
		alncode, alnstdout, alnstderr = run(command)
		
		command = 'eukdetect --mode filter --configfile tests/configfile_for_tests.yml --force'
		filtercode, filterstdout, filterstderr = run(command)
		self.assertTrue(alncode==0 and filtercode == 0, msg=alnstderr + filterstderr)

class _05_cleanup(unittest.TestCase):

	def test_class(self):
		#cleanup
		if os.path.isdir("tests/aln"):
			shutil.rmtree("tests/aln")
		if os.path.isdir("tests/filtering"):
			shutil.rmtree("tests/filtering")
		if os.path.isfile("tests/test_filtered_hits_table.txt"):
			os.remove("tests/test_filtered_hits_table.txt")
		if os.path.isfile("tests/test_filtered_hits_eukfrac.txt"):
			os.remove("tests/test_filtered_hits_eukfrac.txt")
		if os.path.isfile("tests/alignment_commands.txt"):
			os.remove("tests/alignment_commands.txt")
		self.assertTrue(1==1)

if __name__ == '__main__':
	unittest.main()
