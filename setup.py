import setuptools
 
with open("README.md", "r") as fh:
	long_description = fh.read()
 
setuptools.setup(
		name="eukdetect",
		version="v2.0.0",
		author="Abigail Lind",
		description="Detect eukaryotes from shotgun metagenomic data",
		long_description=long_description,
		long_description_content_type="text/markdown",
		license="MIT",
		url="https://github.com/allind/EukDetect",
		entry_points={"console_scripts": [
		"eukdetect = eukdetect.wrappers.main:main",
		"eukdetect-normalize = eukdetect.util.normalize_rpks:main",
	]},
		packages=setuptools.find_packages(),
		package_data={
			"eukdetect.rules": ["*.rules", "*.yaml"],
			"eukdetect.envs":  ["*.yml"],
		},
		include_package_data=True,
		python_requires=">=3.8",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Intended Audience :: Science/Research",
			"Topic :: Scientific/Engineering :: Bio-Informatics",
			"License :: OSI Approved :: MIT License",
			"Programming Language :: Python :: 3",
			"Programming Language :: Python :: 3.8",
			"Programming Language :: Python :: 3.9",
			"Programming Language :: Python :: 3.10",
			"Programming Language :: Python :: 3.11",
			"Programming Language :: Python :: 3.12",
			"Programming Language :: Python :: 3.13",
		],
)
