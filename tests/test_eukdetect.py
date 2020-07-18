#!/usr/bin/env python

import unittest
import os
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
		command = 'eukdetect --mode runall --configfile tests/test_se_configfile.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)

class _03_RunAln(unittest.TestCase):

	def test_class(self):
		command = 'eukdetect --mode aln --configfile tests/test_se_configfile.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)

class _04_RunAlncmd(unittest.TestCase):

	def test_class(self):
		command = 'eukdetect --mode alncmd --configfile tests/test_se_configfile.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)


class _05_RunFilter(unittest.TestCase):

	def test_class(self):
		command = 'eukdetect --mode filter --configfile tests/test_se_configfile.yml --force'
		code, stdout, stderr = run(command)
		self.assertTrue(code==0, msg=stderr)

if __name__ == '__main__':
	unittest.main()
