# Supplementary repository for the FMSI paper --- UNDER UPDATE (Jan 26, 2025)

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

### Experimental results

Our full experimental results are available at [experiments/01_build_and_query_memtime/99_results](experiments/01_build_and_query_memtime/99_results), with one TSV table per dataset (including both original and subsampled reference).

### Reproducing experimental results

We provide a way to run all the experiments reported in the paper at once, using [Snakemake](https://snakemake.readthedocs.io/en/stable/). We also require Rscript to aggregate the results into TSV tables (one per dataset).
After cloning this repository, perform the following preparation steps.

1. Download all (or selected) datasets from [Zenodo record 14722244](https://zenodo.org/records/14722244) into the [data/](data/) directory.
2. Decompress all xz files in [data/](data/).
3. Download and compile the required programs CBL, FMSI, GGCAT, KmerCamel, ProphAsm, SBWT, SSHash, and Wgsim: Go to [experiments/01_build_and_query_memtime/](experiments/01_build_and_query_memtime/) and run all `download_compile_{prog}.sh` scripts. Note that SBWT and SSHash need to be compile separately for $k \ge 32$ or $k > 32$. Also 
CBL requires to be compiled for each value of *k* separately (the range of $k$ for compilation can be limited in [experiments/01_build_and_query_memtime/download_compile_CBL.sh](download_compile_CBL.sh)).
4. The evaluated values of *k*, subsampling rates *r*, and datasets can all be changed in the [Snakefile](experiments/01_build_and_query_memtime/Snakefile) (lines 10-25). Note that running the experiments on all datasets with all combinations of the parameters may take several days of CPU time, a lot of disk space for indexes and FASTA files with SPSS/MS (about 1 TB), and some computations may require 200 GB or more. The number of cores provided to Snakemake can be changed in the [Makefile](experiments/01_build_and_query_memtime/Makefile) (default = 2). Note also that running SBWT requires substantial disk space for temporary files, especially on the *E. coli* pan-genomes, HG T2T assembly, and HG Illumina reads  (e.g., on the whole *E. coli* pan-genome, it used 84 GB of disk space in our benchmarks). 

Then, the experiments are run with `make`. Use `make test` to run it only on *S. pneumoniae* pan-genome with $k=31$.

## Figures

TODO: UPDATE: The [figures/](figures/) directory contains Fig. 2 and the associated supplementary plots, including the script used for their plotting.

## Contact

* Ondřej Sladký (ondra.sladky@gmail.com)
* Pavel Veselý (vesely@iuuk.mff.cuni.cz)
* Karel Břinda (karel.brinda@inria.fr)
