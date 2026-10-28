#!/usr/bin/env python
"""Import ete3 on Python versions that no longer ship the `cgi` module.

ete3 3.1.3 imports its web plugin from package __init__, and that plugin does
`import cgi`. PEP 594 removed `cgi` from the standard library in Python 3.13,
so `from ete3 import NCBITaxa` raises ModuleNotFoundError there even though
nothing in EukDetect touches the web plugin.

The plugin uses exactly one name from the module, `cgi.parse_qs`, which has
always been an alias for `urllib.parse.parse_qs`. Supplying that under the
name `cgi` before ete3 is imported satisfies the import with stdlib code and
removes the need for an external backport package.

Import NCBITaxa from here rather than from ete3 directly.
"""
import importlib
import sys
import types
from urllib.parse import parse_qs, parse_qsl


def install_cgi_shim():
	"""Register a minimal stand-in for the removed `cgi` module.

	Returns True if a shim was installed, False if a real `cgi` was importable.
	"""
	if "cgi" in sys.modules:
		return False
	try:
		importlib.import_module("cgi")
		return False
	except ModuleNotFoundError:
		pass

	shim = types.ModuleType("cgi")
	shim.__doc__ = (
		"Minimal stand-in for the removed standard library cgi module, "
		"installed by eukdetect.util.ete3_compat so that ete3 can be imported "
		"on Python 3.13 and later. Only parse_qs and parse_qsl are provided."
	)
	shim.parse_qs = parse_qs
	shim.parse_qsl = parse_qsl
	sys.modules["cgi"] = shim
	return True


install_cgi_shim()

from ete3 import NCBITaxa  # noqa: E402

__all__ = ["NCBITaxa", "install_cgi_shim"]
