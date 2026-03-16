#  EukDetect
 
 EukDetect is a bioinformatics tool for detecting eukaryotic organisms in metagenomic sequencing data. It uses a curated database of marker genes to identify eukaryotic species present in microbial communities.
 
## Features
 
 - Detect eukaryotic organisms from shotgun metagenomic sequencing data
 - Support for both paired-end and single-end sequencing data
 - Process individual samples or batches in parallel
 - Estimates absolute and relative abundance of detected eukaryotes
 - Comprehensive quality filtering and marker gene coverage analysis
 
## Installation
 
### Option 1: Install from Bioconda (Recommended)
 
 TBD
 
### Option 2: Install from GitHub
 
 **Installation steps:**
 
 1. Clone the repository:
 ```
 git clone https://github.com/allind/EukDetect.git
 cd EukDetect
 ```
 
 2. Create the conda environment:
 ```
 conda env create -f eukdetect/envs/eukdetect2_environment.yml
 conda activate eukdetect
 ```
 
 3. Install EukDetect in development mode:
 ```
 pip install -e .
 ```
 
 4. Verify installation:
 ```
 eukdetect --help
 ```
 
### Database Installation
 
 **Download EukDetect database from Figshare**
 
 Download and unpack the EukDetect database (eukdetect2_database.tar.gz) from the [Figshare repository](https://doi.org/10.6084/m9.figshare.31360168).
 
 The previous versions of the EukDetect database are [still available](https://doi.org/10.6084/m9.figshare.12670856.v4), but are not compatible with Eukdetect2. #update this
 
 ```
 mkdir eukdb
 cd eukdb
 wget [tbd]
 tar -zxvf [tbd]
 rm [tbd]
 ```
 
 The uncompressed database folder is around 7 GB.
 
## Quick Start
 
### Single Sample
 
 Process a single paired-end sample:
 ```
 eukdetect single \
   -1 sample_R1.fastq.gz \
   -2 sample_R2.fastq.gz \
   -n sample_name \
   --outdir results/ \
   --database /path/to/eukdetect/database \
   --database-prefix eukdb \
   --cores 16
 ```
 
 Process a single-end sample:
 ```
 eukdetect single \
   -1 sample.fastq.gz \
   -n sample_name \
   --outdir results/ \
   --database /path/to/eukdetect/database \
   --database-prefix eukdb \
   --cores 16
 ```
 
 Single-sample cores are used to run multi-threaded bowtie2 and are not used for other parts of the pipeline.
 
 There is currently no support for a mixture of paired and single end reads for a sample. 
 
### Multiple Samples (Batch local mode)
 
 1. Create a tab-separated samples file (`samples.tsv`):
 
 For paired-end data:
 ```
 sample_name	reads1	reads2
 sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
 sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz
 sample3	/path/to/sample3_R1.fastq.gz	/path/to/sample3_R2.fastq.gz
 ```
 
 For single-end data:
 ```
 sample_name	reads1
 sample1	/path/to/sample1.fastq.gz
 sample2	/path/to/sample2.fastq.gz
 sample3	/path/to/sample3.fastq.gz
 ```
 
 2. Run batch mode:
 ```
 eukdetect batch \
   --samples samples.tsv \
   --outdir results/ \
   --database /path/to/eukdetect/database \
   --database-prefix eukdb \
   --cores 10
 ```
 
 This will process 10 samples in parallel on your local machine. Bowtie2 will use 1 core by default.
 
## Usage
 
### Command-line Interface
 
 EukDetect provides two main commands:
 
#### Single Mode
 Process individual samples. Use this for
 - Running one sample at a time
 - HPC job arrays (one job per sample)
 
 ```
 eukdetect single [options]
 ```
 
 **Required arguments:**
 - `-1, --reads1`: Forward reads (R1) or single-end reads (required)
 - `--outdir, -o`: Output directory (required)
 - `--database, -d`: Path to EukDetect database directory (required)
 
 **Optional arguments:**
 - `-2, --reads2`: Reverse reads (R2) for paired-end data
 - `-n, --sample-name`: Sample name (auto-detected from filename if not provided)
 - `--database-prefix`: Database file prefix (default: eukdb)
 - `--cores, -c`: Number of CPU threads for alignment (default: 1)
 - `--readlen`: Read length (auto-detected if not provided)
 - `--mode`: Analysis mode (default: all)
 - `--dry-run`: Preview commands without executing
 
#### Batch Mode
 Process multiple samples in parallel locally.
 
 ```
 eukdetect batch [options]
 ```
 
 **Required arguments:**
 - `--samples`: Tab-separated file with sample information (required)
 - `--outdir, -o`: Output directory (required)
 - `--database, -d`: Path to EukDetect database directory (required)
 
 **Optional arguments:**
 - `--database-prefix`: Database file prefix (default: eukdb)
 - `--cores, -c`: Number of samples to process in parallel (default: 1)
 - `--readlen`: Read length in bp (auto-detected if not provided)
 - `--mode`: Analysis mode (default: all)
 - `--dry-run`: Preview workflow without executing
 
 Unlike single batch mode, cores are used to parallelize samples and the alignment step uses a single core per sample.
 
### Analysis Modes
 
 EukDetect supports four analysis modes via the `--mode` option:
 
 - `all` (default): Run complete pipeline (alignment + filtering + analysis)
 - `aln`: Run alignment step only
 - `analyze`: Run filtering and analysis only (requires existing alignments)
 - `printaln`: Generate file with alignment commands for manual execution
 
 **Examples:**
 
 ```
# Run only alignment
 eukdetect single -1 R1.fq.gz -2 R2.fq.gz -n sample -o out/ -d db/ --mode aln --cores 16
 
# Run analysis on existing alignments
 eukdetect single -1 R1.fq.gz -2 R2.fq.gz -n sample -o out/ -d db/ --mode analyze
 
# Generate alignment commands for manual execution
 eukdetect single -1 R1.fq.gz -2 R2.fq.gz -n sample -o out/ -d db/ --mode printaln
 ```
 
## Output Files
 
 EukDetect creates the following output directory structure:
 
 ```
 results/
 ├── configs/
 │   └── config_sample1.yml          
 ├── logs/
 │   └── snakemake_sample1_*.log      
 ├── aln/
 │   └── sample1_aln_q10_lenfilter.sorted.bam 
 ├── filtering/
 │   ├── sample1_aln_q10_lenfilter.dedup.bam
 │   └── sample1_read_counts_and_mismatches.txt
 │   └── sample1_all_hits_table.txt
 ├── sample1_filtered_hits_table.txt
 └── sample1_filtered_hits_eukfrac.txt
 ```
 
 **`{sample}_filtered_hits_table.txt`** reports for each detected taxon:
 - Taxonomic name, rank, lineage, and NCBI taxonomy ID
 - Number of marker genes with aligned reads
 - Total number of reads aligning to marker genes
 - Percent_observed_markers: percentage of marker genes detected
 - Total_marker_coverage: percentage of bases covered in observed markers
 - Percent_identity: average percent identity across aligned reads
 
 **`{sample}_filtered_hits_eukfrac.txt`** reports:
 - RPKS (Reads Per Kilobase of Sequence): absolute abundance metric for **species-level taxa only**
 - Relative_abundance (EukFrac): relative abundance compared to other eukaryotes at all taxonomic levels
 - Total reads: number of reads aligning to markers at each taxonomic level
 
##  Important Considerations
 **Eukfrac interpretation:** The EukFrac (relative abundance) metric is relative only to other eukaryotes, not to bacteria or archaea. Always consider RPKS alongside EukFrac when interpreting results.

 **Filtering threshholds:** By default, EukDetect removes taxa with fewer than 4 reads aligning to fewer than 2 marker genes. Unfiltered results are available in the `filtering/` directory.
 
 **Read length:** EukDetect supports reads over 75 base pairs long.
 
 **Paired-end data:** Forward and reverse read files must have the same number of reads (properly paired).
 
 **Mixed data:** Cannot mix single-end and paired-end samples in the same batch run. Run them separately.
 
 **File naming:** For batch mode, all samples must use consistent file extensions and naming patterns.

 **Input file directory:** All input FASTQ files for a single batch run must reside in the same directory. EukDetect infers one shared `fq_dir` from the first sample and constructs all other file paths as `{fq_dir}/{sample_name}{suffix}`. Samples stored in different directories cannot be processed in a single batch run — symlink them into a common directory first, or run them as separate `eukdetect single` invocations.

**Comparing across samples:** RPKS cannot be directly compared between samples without first normalizing by library size.
  
## Normalizing RPKS Across Samples

RPKS values from EukDetect are not normalized by total sequencing depth. This means RPKS values cannot be directly compared across samples with different library sizes. To enable cross-sample comparisons and to combine RPKS and eukfrac data across samples into one file, use the provided `eukdetect-normalize` command

The normalization script produces a combined table with additional columns:

- **Sample**: Sample name
- **Library_size**: Total reads in the library
- **RPSKM**: Normalized abundance (RPKS per million total reads)

RPKSM values are only calculated for species-level taxa. Higher taxonomic levels show "NA" for RPKSM.

**Usage:**
```
eukdetect-normalize \
  --eukfrac results/*_filtered_hits_eukfrac.txt \
  --library-sizes library_sizes.tsv \
  --output all_samples_normalized.tsv
```

**Preparing Library Sizes File**
Create a tab-separated file with sample names and total read counts:

```
sample_name	total_reads
sample1	50000000
sample2	75000000
sample3	100000000
```

**To count reads from FASTQ files:**

```
# Create header
echo -e "sample_name\ttotal_reads" > library_sizes.tsv

# Count reads for each sample
for sample in sample1 sample2 sample3; do
    count=$(zcat ${sample}_1.fastq.gz | echo $((`wc -l`/4)))
    echo -e "${sample}\t${count}" >> library_sizes.tsv
done
```


 
 
## Testing
 
 Run the test suite to verify your installation:
 
 ```
 pytest tests/ -v
 ```
 
## Citation
 
 If you use EukDetect, please cite the [paper in _Microbiome_](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01015-y):
 
 Lind, A.L., Pollard, K.S. Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing. Microbiome 9, 58 (2021).
 
## Taxonomy Database Issues
 
 The EukDetect pipeline uses the ete3 package to interface with the NCBI taxonomy database. The database uses the NCBI taxonomy release from early 2026. The Figshare repository includes both the taxdump file and the ete3 sqlite database.
 
 If you encounter errors from the ete3 package, you may need to regenerate the taxonomy database:
 
 ```
 conda activate eukdetect
 ```
 
 Open a Python console and run:
 
 ```
 from ete3 import NCBITaxa
 ncbi = NCBITaxa()
 ncbi.update_taxonomy_database(taxdump_file="taxdump_1_14_2020.tar.gz")
 exit()
 ```
 
 Move the newly created `taxa.sqlite` and `taxa.sqlite.traverse.pkl` files from `~/.etetoolkit/` to the EukDetect database folder.
 
## License
 
 EukDetect is distributed under the MIT License. See LICENSE file for details.
