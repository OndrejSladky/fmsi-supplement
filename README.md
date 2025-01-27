# Supplementary repository for the FMSI paper

Ondřej Sladký, Pavel Veselý, Karel Břinda:
**"FroM Superstring to Indexing: a space-efficient index for unconstrained *k*-mer sets using the Masked Burrows-Wheeler Transform (MBWT)"**, [preprint at *bioRxiv*](https://www.biorxiv.org/content/10.1101/2024.10.30.621029), 2024.

### Citation

```
@article {Sladky2024.10.30.621029,
	author = {Sladk{\'y}, Ond{\v r}ej and Vesel{\'y}, Pavel and B{\v r}inda, Karel},
	title = {FroM Superstring to Indexing: a space-efficient index for unconstrained k-mer sets using the Masked Burrows-Wheeler Transform (MBWT)},
	elocation-id = {2024.10.30.621029},
	year = {2024},
	doi = {10.1101/2024.10.30.621029},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/11/03/2024.10.30.621029},
	eprint = {https://www.biorxiv.org/content/early/2024/11/03/2024.10.30.621029.full.pdf},
	journal = {bioRxiv}
}

```

## Table of Contents

<!-- vim-markdown-toc GFM -->

* [Experimental evaluation](#experimental-evaluation)
* [Figures](#figures)

<!-- vim-markdown-toc -->


## Experimental evaluation

### Benchmark datasets

* Two *E. coli* pan-genomes, obtained as the union of the E. coli genomes from the 661k collection. One contains all genomes (without quality filtering) and for the other (HQ) we applied high-quality filtering.
* *S. pneumoniae* pan-genome: 616 genomes, as provided in RASE DB S. pneumoniae https://github.com/c2-d2/rase-db-spneumoniae-sparc/
* *SARS-CoV-2* pan-genome, downloaded from GISAID https://gisaid.org/ (access upon registration) on Jan 25, 2023 (GISAID version 2023/01/23, 14,682,066 genomes, 430 Gbp).
* Metagenomic sample SRS063932 (Illumina raw reads) of human microbiome with accession SRX023459, download from https://www.hmpdacc.org/hmp/HMASM/. The fastq files were converted to FASTA files using `seqtk seq -A -C`.
* Human RNA-seq Illumina raw reads with accession SRX348811, downloaded using the prefetch tool from the SRA toolkit and then converted into the FASTA format by `fastq-dump --split-3 --fasta`.
* Human genome Illumina raw reads with accession SRX016231, downloaded using the prefetch tool from the SRA toolkit and then converted into the FASTA format by `fastq-dump --split-3 --fasta`.
* Human genome assembly chm13v2.0 (T2T), downloaded from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.
* Two MiniKraken datasets (4GB and 8GB), downloaded from https://ccb.jhu.edu/software/kraken/, with the 31-mers dumped using Jellyfish 1.1.12.

We subsampled some of the datasets (pan-genomes, metagenomic sample, and HG T2T assembly) by choosing 10% randomly chosen distinct canonical $k$-mers. 

All of the datasets are provided at [Zenodo record 14722244](https://zenodo.org/records/14722244).

For each of the pan-genomes, we add one reference genome to the [data/](data/) directory, downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/). Other datasets also have a "reference genome" available in [data/](data/) (for HG T2T assembly and HG Illumina reads, it is the HG T2T assembly). These reference genomes are used to generate positive streaming queries using [Wgsim](https://github.com/lh3/wgsim).

For generating negative membership queries to these datasets, isolated queries, we used chromosome 1A of \emph{T. aestivum} assembly `GCF\_018294505.1`, downloaded from [NCBI](https://www.ncbi.nlm.nih.gov), provided at [data/Triticum_aestivum.IWGSC.dna.chromosome.1A.fa.xz](data/Triticum_aestivum.IWGSC.dna.chromosome.1A.fa.xz).

### Experimental results on membership queries

Our full experimental results for membership queries are available at [experiments/01_build_and_query_memtime/99_results](experiments/01_build_and_query_memtime/99_results), with one TSV table per dataset (including both original and subsampled reference).

### Reproducing experimental results on membership queries

We provide a way to run all the experiments with membership queries reported in the paper, using [Snakemake](https://snakemake.readthedocs.io/en/stable/). We evaluate CBL, FMSI, SBWT, and SSHash. We also require [seqtk](https://github.com/lh3/seqtk) to process FASTA files and Rscript to aggregate the results into TSV tables (one per dataset).
After cloning this repository, perform the following preparation steps.

1. Download all (or selected) datasets from [Zenodo record 14722244](https://zenodo.org/records/14722244) into the [data/](data/) directory.
2. Decompress all xz files in [data/](data/).
3. Download and compile the required programs CBL, FMSI, GGCAT, KmerCamel, ProphAsm, SBWT, SSHash, and Wgsim: Go to [experiments/01_build_and_query_memtime/](experiments/01_build_and_query_memtime/) and run all `download_compile_{prog}.sh` scripts. Note that SBWT and SSHash need to be compile separately for $k \ge 32$ or $k > 32$. Also 
CBL requires to be compiled for each value of *k* separately (the range of $k$ for compilation can be limited in [experiments/01_build_and_query_memtime/download_compile_CBL.sh](download_compile_CBL.sh)).
4. The evaluated values of *k*, subsampling rates *r*, and datasets can all be changed in the [Snakefile](experiments/01_build_and_query_memtime/Snakefile) (lines 10-25). 

Then, the experiments are run with `make` in [experiments/01_build_and_query_memtime/](experiments/01_build_and_query_memtime/). Use `make test` to run it only on *S. pneumoniae* pan-genome with $k=31$.

Notes:
- The number of CPU cores provided to Snakemake can be changed in the [Makefile](experiments/01_build_and_query_memtime/Makefile) (default = 2).
- Running the experiments on all datasets with all combinations of the parameters may take several days of CPU time, about 1 TB of disk space for indexes and FASTA files with SPSS/MS, and some computations may require 200 GB or more.  
- Running SBWT requires substantial disk space for temporary files, especially on the *E. coli* pan-genomes, HG T2T assembly, and HG Illumina reads  (e.g., on the whole *E. coli* pan-genome, it used 84 GB of disk space in our benchmarks); hence, the Snakemake resource parameter `sbwtlimit` (provided as a command-line parameter in Makefile) limit the number of parallel executions of `sbwt build` (default = 1).

### Experiments with dictionary queries

We also provide a similar experimental evaluation for dictionary queries, namely the lookup queries as called in the [SSHash paper](https://doi.org/10.1093/bioinformatics/btac245), which evaluate a minimal perfect hash function (MPHF) for the indexed k*-mers. We only compared FMSI and SSHash as other tools do not provide dictionary functionality. Our results for the MiniKraken datasets are provided in [experiments/03_lookup_queries/99_results](experiments/03_lookup_queries/99_results).

The setup to run the experiments is very similar to the setup for membership queries; in fact, we reuse some indexes from the previous experiments. The FMSI index need to be computed again, on the masked superstring with the minimum number of ones. SSHash index construction and queries are the same as for membership queries.

The experiments for lookup queries are run with `make` [experiments/03_lookup_queries/](experiments/03_lookup_queries/). Use `make test` to run it only on *S. pneumoniae* pan-genome with $k=31$.

## Figures

The [figures/](figures/) directory contains source files for Figures 2-4, Table 2, and supplementary Table 3, including the script used for their plotting.

## Contact

* Ondřej Sladký (ondra.sladky@gmail.com)
* Pavel Veselý (vesely@iuuk.mff.cuni.cz)
* Karel Břinda (karel.brinda@inria.fr)
