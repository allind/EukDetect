# EukDetect Secondary Hits

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

Download and unpack the EukDetect database (eukdetect_database_v1.tar.gz) from the [Figshare repository](https://doi.org/10.6084/m9.figshare.12670856).

```
wget https://ndownloader.figshare.com/files/24012515 
tar -zxvf 24012515
rm 24012515
```

The uncompressed database folder is 2.6 Gb in size.

**Create conda environment and install EukDetect**

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

A schematic of the eukdetect pipeline and the files created in the pipeline can be found in [eukdetect_pipeline_schematic.pdf](https://github.com/allind/EukDetect/blob/master/eukdetect_pipeline_schematic.pdf).

There are *four* eukdetect modes, invoked with `eukdetect --mode`. All modes require a eukdetect config file as described above.

The `runall` mode runs the entire pipeline. 

The `aln` mode runs just the bowtie2 alignment step. 

The `filter` mode runs everything downstream of the alignment. The filter mode can only be run if the alignment step has been completed.

The `alncmd` mode  creates a file in the output directory specified in the configfile called alignment_commands.txt. These commands can be run on a compute cluster either sequentially or as a job array. The alignment step of eukdetect is the most computationally intensive step of eukdetect, and this mode is intended for users to run the alignments on a high performance compute cluster if desired.

Examples of eukdetect usage:

```
eukdetect --mode runall --configfile [config file] --cores [cores]
eukdetect --mode aln --configfile [config file] --cores [cores]
eukdetect --mode filter --configfile [config file] --cores [cores]
eukdetect --mode alncmd --configfile [config file] --cores [cores]
```

**Important info**

Currently, EukDetect only supports analysis of reads that are **over** 75 base pairs long.

If you are running EukDetect in paired-end mode, the number of reads in both files has to match. If there are different numbers of reads in the forward and reverse files, EukDetect will fail at the alignment phase.

A list of currently detectable EukDetect species is in Table S2 in the supplementary material of the [biorxiv manuscript](https://www.biorxiv.org/content/10.1101/2020.07.22.216580v1).


## Output file descriptions

A schematic of the eukdetect pipeline and the files created in the pipeline can be found in [eukdetect_pipeline_schematic.pdf](https://github.com/allind/EukDetect/blob/master/eukdetect_pipeline_schematic.pdf).

The main output of EukDetect are the files `{output_directory}/{samplename}_stats_per_filtered_taxid.txt` and `{output_directory}/{samplename}_hit_taxonomy_filterpass.txt`. 

`{samplename}_stats_per_filtered_taxid.txt` reports: 
  * the observed taxon's name
  * how many marker genes had >1 aligned read
  * the overall number of reads aligning to marker genes from this taxon
  
It also reports three statistics calculated from these numbers, which are:
  * Percent_observed_markers: the percentage of the total marker genes of that species that are observed
  * Total_marker_coverage: of the markers that are observed for that taxon, what percentage of the bases in that marker gene have one or more aligned reads
  * Percent_identity: the percent identity calculated across all reads aligning to that taxon's marker gene

`{samplename}_hit_taxonomy_filterpass.txt` is a taxonomy tree of all observed taxa, and it reports:
  * Markers_Obs: The number of markers observed at that taxonomic node (includes all markers in node children)
  * Total_Markers: The total number of markers that can be observed (includes all markers in node children)
  * Percent_Makers_Obs: Percentage of total_markers observed at this node
  * Percent_ID: percent identity shown in `{samplename}_stats_per_filtered_taxid.txt`. Only calculated if reads are assigned to this exact node, not to node children.
  * Marker_read_count: number of reads aligning to markers at this node (including markers of node children).
  * Rank: NCBI taxonomy rank of node

EukDetect removes all taxa that have fewer than 4 reads that align to fewer than 2 marker genes. This is the minimum amount of evidence we recommend for determining if a eukaryotic species is present. However, the same information as the main output files without any filtering is located in the `{output_directory}/filtering/` folder. Information about reads aligning to each marker is in `{output_directory}/filtering/{samplename}_read_counts_and_mismatches.txt`.

Alignments to the databse that have been length and quality filtered are in a coordinate-sorted bam file in the `{output_directory}/aln/` folder. Alignments that additionally have low-complexity and duplicate reads removed are in a bam file in `{output_directory/filtering/`.

## Citation

If you use EukDetect, please cite our [preprint on biorxiv](https://www.biorxiv.org/content/10.1101/2020.07.22.216580v1).

<a id="1"></a> 
Lind, A. L. and K.S. Pollard. (2020).
Accurate and sensitive detection of microbial eukaryotes from metagenomic shotgun sequencing data.
bioRxiv. 2020.07.22.216580. doi: https://doi.org/10.1101/2020.07.22.216580


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
