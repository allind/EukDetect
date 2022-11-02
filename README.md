# EukDetect

## Installation

**Install conda**

If you do not have the conda package manager installed already, follow the [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to install Miniconda.

**Download this repository**

Download this repository & download the Fishare repo
```
git clone https://github.com/allind/EukDetect.git
cd EukDetect
```

**Download EukDetect database from Figshare**

***Update 4/22/2022 - New database uploaded on Figshare with files required for estimating eukaryotic relative abundance.***

***NEW DATABASE 1/23/2021 - MUCH LARGER TAXONOMIC COVERAGE***

Download and unpack the EukDetect database (eukdetect_database_v2.tar.gz) from the [Figshare repository](https://doi.org/10.6084/m9.figshare.12670856.v8)

The previous version of the EukDetect database (from NCBI genomes only, no chloroplastids and no metazoans other than worms) is [still available](https://doi.org/10.6084/m9.figshare.12670856.v4).

```
mkdir eukdb
cd eukdb
wget https://ndownloader.figshare.com/files/34885596
tar -zxvf 34885596
rm 34885596
```

The uncompressed database folder is 2.6 Gb in size.

**Create conda environment and install EukDetect**

Run this code inside the EukDetect main folder, where the environment.yml file is located.

```
conda env update --name eukdetect -f environment.yml
conda activate eukdetect
# install eukdetect
python setup.py install
```

**Test installation**

To test your installation, edit the file `configfile_for_tests.yml` with the path to the installation directory and the path to the EukDetect database.

From within the EukDetect installation directory, run:

```
python tests/test_eukdetect.py
```

## Usage

**Edit the config file**

Copy the `default_configfile.yml` to `your_configfile.yml`. Change all parameters in the config file as described.

If you don't know what the length of your reads is, this is a handy one-liner to estimate it: `gzip -dc {file.fastq.gz} | head -n 10000 | awk '{ if (NR%4==2){count++; bases += length}} END{printf "%3.0f\n", bases/count}'`

**Run modes**

A schematic of the eukdetect pipeline and the files created in the pipeline can be found in [the pipeline schematic pdf](https://github.com/allind/EukDetect/blob/master/EukDetect_pipeline_schematic.pdf).

There are *four* eukdetect modes, invoked with `eukdetect --mode`. All modes require a eukdetect config file as described above.

The `runall` mode runs the entire pipeline. 

The `aln` mode runs just the bowtie2 alignment step. 

The `filter` mode runs everything downstream of the alignment. The filter mode can only be run if the alignment step has been completed.

The `printaln` mode  creates a file in the output directory specified in the configfile called alignment_commands.txt. These commands can be run on a compute cluster either sequentially or as a job array. The alignment step of eukdetect is the most computationally intensive step of eukdetect, and this mode is intended for users to run the alignments on a compute cluster if desired.

Examples of eukdetect usage:

```
eukdetect --mode runall --configfile [config file] --cores [cores]
eukdetect --mode aln --configfile [config file] --cores [cores]
eukdetect --mode filter --configfile [config file] --cores [cores]
eukdetect --mode printaln --configfile [config file] --cores [cores]
```


**Running snakemake directly**

EukDetect can also be run directly as a snakemake workflow using the `rules/eukdetect.rules` file, specifying either `runall`, `printaln`, `aln`, or `filter` as the target rule. If you routinely run snakemake jobs on a cluster and wish to run the entire EukDetect pipeline on it, this is the recommended option. If you're running into issues with the eukdetect python package this is also the recommended running option. Running snakemake directly means there are fewer checks to make sure the input and output are correct.

Examples:
```
snakemake --snakefile rules/eukdetect.rules --configfile [config file] --cores [cores] runall
```

**Important info**

Currently, EukDetect only supports analysis of reads that are **over** 75 base pairs long.

If you are running EukDetect in paired-end mode, the number of reads in both files has to match. If there are different numbers of reads in the forward and reverse files, EukDetect will fail at the alignment phase.

A list of currently detectable EukDetect species is in Table S2 in the supplementary material of the [biorxiv manuscript](https://www.biorxiv.org/content/10.1101/2020.07.22.216580v1).


## Output file descriptions

A schematic of the eukdetect pipeline and the files created in the pipeline can be found in [eukdetect_pipeline_schematic.pdf](https://github.com/allind/EukDetect/blob/master/eukdetect_pipeline_schematic.pdf).

The main output of EukDetect are the files `{output_directory}/{samplename_filtered_hits_table.txt` and `{output_directory}/{samplename}_filtered_hits_eukfrac`. 

`{samplename}_filtered_hits_table.txt` reports: 
  * the observed taxon's name, taxonomic rank, taxonomic lineage, and NCBI taxonomy ID (all from the NCBI taxonomy database)
  * how many marker genes had >1 aligned read
  * the overall number of reads aligning to marker genes from this taxon
  
It also reports three statistics calculated from these numbers, which are:
  * Percent_observed_markers: the percentage of the total marker genes of that species that are observed
  * Total_marker_coverage: of the markers that are observed for that taxon, what percentage of the bases in that marker gene have one or more aligned reads
  * Percent_identity: the percent identity calculated across all reads aligning to that taxon's marker gene

`{samplename}_filtered_hits_eukfrac.txt` reports the taxonomic lineage of all hits, along with a calculation of ```RPKS (Reads Per Kilobase of Sequence```, which is related to the absolute abundance of taxa, and ```EukFrac (Eukaryotic Fraction)```, which is related to the abundance of taxa relative to other eukaryotes. RPKS is calculated by dividing the reads aligned to markers by the length in Kb of all markers for a species. The EukFrac relative abundance is calculated for each of 6 taxonomic levels - phylum, order, class, family, genus, and species. If any hit does not have one of these taxonomic levels (which occasionally happens in the NCBI taxonomy database), the EukFrac is not calculated for that taxonomic level.

**Important Note**:
It is **very important** to note that the EukFrac metric is only relative to **other eukaryotes** and not to bacteria or archaea, which almost always make up the majority of a microbiome sequencing library. The EukFrac metric must be considered alongside a proxy for the absolute abundance of taxa, which is RPKS. I strongly discourage directly comparing changes in EukFrac alone between samples (such as is commonly done with stacked bar charts) because the fraction of the sequencing library comprised by eukaryotes is so small and these measures are noisy.

EukDetect removes all taxa that have fewer than 4 reads that align to fewer than 2 marker genes. This is the minimum amount of evidence we recommend for determining if a eukaryotic species is present. However, the same information as the main output files without any filtering is located in the `{output_directory}/filtering/` folder. Information about reads aligning to each marker is in `{output_directory}/filtering/{samplename}_read_counts_and_mismatches.txt`.

Alignments to the database that have been length and quality filtered are in a coordinate-sorted bam file in the `{output_directory}/aln/` folder. Alignments that additionally have low-complexity and duplicate reads removed are in a bam file in `{output_directory/filtering/`.

## Citation

If you use EukDetect, please cite the [paper in _Microbiome_](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01015-y).

<a id="1"></a> 
Lind, A.L., Pollard, K.S. Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing. Microbiome 9, 58 (2021).


## Taxonomy database version

**This section is only necessary to read if you are getting errors from the ete3 package.**

The EukDetect pipeline uses the ete3 package to interface with the NCBI taxonomy database. In order for the EukDetect pipeline to run correctly, ete3 must use the NBCI taxonomy database corresponding to the January 14, 2020 taxonomy release. The Figshare repository for the EukDetect database has both the taxdump_1_14_2020.tar.gz file, as well as the ete3 sqlite database (in files taxa.sqlite and taxa.sqlite.traverse.pkl). EukDetect uses the database located in this folder.

It may be possible that updates to the ete3 package, or differences in operating systems, might result in ete3 not being able to correctly parse the database. If this happens, users will need to create their own taxonomy database from taxdump_1_14_2020.tar.gz. By default, ete3 saves the taxonomy database in ~/.etetoolkit in the home directory of the user, so this file will need to be moved into the EukDetect database directory. If you use ete3 outside of EukDetect, you may run the ete3 function update_taxonomy_database(), which will download the newest NCBI taxonomy database and overwrite the taxonomy database in ~/.etetoolkit. Therefore, this file needs to be moved to a location where it will not be overwritten.

**Solution**

If you use the ete3 package outside of EukDetect and you have a specific version of the NCBI taxonomy database installed, copy ~/.etetoolkit/taxa.sqlite and ~/.etetoolkit/taxa.sqlite.traverse.pkl to temporary files that you will restore afterwards. Otherwise, this process will overwrite any existing files and you will lose the data.

```
conda activate eukdetect
```

open a Python console and run the following code:

```
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database(taxdump_file="taxdump_1_14_2020.tar.gz")
exit()
```

Now, remove the taxa.sqlite and taxa.sqlite.traverse.pkl file from the database folder, and find the newly created taxa.sqlite and taxa.sqlite.traverse.pkl files. This file will be located in your home directory in ~/.etetoolkit/. Move these file into the EukDetect database folder. If you moved your existing taxa.sqlite files to temporary files, restore them to ~/.etetoolkit.
