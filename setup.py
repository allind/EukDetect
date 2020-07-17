import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

with open('LICENSE') as f:
	license=f.read()

setuptools.setup(
		name="EukDetect",
		version="0.0.1",
		author="Abigail Lind",
		author_email="allind89@gmail.com",
		description="Detect eukaryotes from shotgun metagenomic data",
		long_description=long_description,
		long_description_content_type="text/markdown",
		license=license,
		url="https://github.com/allind/EukDetect.git",
		entry_points={"console_scripts": ["eukdetect = eukdetect.runall:main"]},
		packages=setuptools.find_packages()
)

