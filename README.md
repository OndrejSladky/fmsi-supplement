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

* *E. coli* pan-genome, obtained as the union of the genomes from the [661k collection](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001421), downloaded from [Phylogenetically compressed 661k collection](https://zenodo.org/records/4602622)
  - Available at [Zenodo](https://zenodo.org/records/13997398)
* *S. pneumoniae* pan-genome - 616 genomes, as provided in [RASE DB *S.
  pneumoniae*](https://github.com/c2-d2/rase-db-spneumoniae-sparc/)
  - *k*-mers were collected and stored in the form of simplitigs (ProphAsm
    v0.1.1, k=32, NS: 158,567, CL: 14,710,895 bp, #kmers: 9,795,318 32-mers)
  - The resulting file:
    [data/spneumo_pangenome_k32.fa.xz](data/spneumo_pangenome_k32.fa.xz)
* *SARS-CoV-2* pan-genome - downloaded from [GISAID](https://gisaid.org/)
  (access upon registration) on Jan 25, 2023 (GISAID version 2023_01_23,
  14,682,066 genomes, 430 Gbp)
  - *k*-mers were collected using JellyFish 2 (v2.2.10, 11,701,570 32-mers) and
    stored in the form of simplitigs (ProphAsm v0.1.1, k=32, NS: 345,866, CL:
    22,423,416 bp, #kmers: 11,701,570 32-mers)
  - The resulting file:
    [data/sars-cov-2_pangenome_k32.fa.xz](data/sars-cov-2_pangenome_k32.fa.xz)

For each of the pan-genomes, we add one reference genome to the [data/](data/) directory, downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/). These reference genomes are used to generate positive streaming queries using [Wgsim](https://github.com/lh3/wgsim).

For generating negative membership queries to these datasets, we used a 2MB prefix of the FASTA file for chromosome 1 of *H. sapiens* genome (`GRCh38.p14 Primary Assembly`, `NC_000001.11`), downloaded from [NCBI](https://www.ncbi.nlm.nih.gov); see  [data/GRCh38.p14.chromosome1.prefix2M.fasta.xz](data/GRCh38.p14.chromosome1.prefix2M.fasta.xz)


### Reproducing experimental results

After cloning this repository, run the following to download all the dependencies.

```bash 
git submodule update --init
```
After that, CBL, SBWT, BWA, FMSI, KmerCamel, ProphAsm, and Wgsim (the submodules) need to be compiled, as described in each of these repositories.
We note that CBL need to be compiled for each value of *k* separately, and we provide script [`compileCBL.sh`](compileCBL.sh) which compiles CBL for $k = 15, 23,$ and 31 with appropriate parameters.

To run the experiments on membership queries for one conmbination of dataset, value of *k*, and subsampling rate, execute
```bash
cd experiments/
./run_experiment_memtime.sh <prefix> <k> <rate>
```
where `<prefix>` is the dataset, e.g., `<prefix> = spneumo_pangenome_k32` or `sars-cov-2_pangenome_k32` or `escherichia_coli.k63` (the script assumes that files `../data/<prefix>.fa.xz` and `../data/<prefix>-refGenome.fa` exist),
`<k>` is the *k*-mer size, and `<rate>` is the rate to subsample *k*-mers (e.g., rate could be 0.1 for sampling 10% distinct *k*-mers; use 1.0 for no subsampling). For example, `./run_experiment_memtime.sh spneumo_pangenome_k32 31 0.1`.

The script creates log files in the `experiments/` directory, each containing the results of /usr/bin/time in a tsv format (with a header line specifying the command that was run).

**Warning:** running SBWT requires substantial disk space for temporary files, especially on the *E. coli* dataset, where it used 84 GB of disk space in our benchmarks.

We also provide a way to run the experiments reported in the paper at once (i.e., run it with multiple datasets, values of *k*, and subsampling rates), using [Snakemake](https://snakemake.readthedocs.io/en/stable/). We also require Rscript to aggregate the results into a tsv table.

First, the datasets evaluated need to be subsampled using script `run_subsampling.sh`, which gets dataset name (without extension .fa.xz) as a parameter. One can specify desired subsampling rates and values of *k* inside `run_subsampling.sh`. This creates compressed FASTA files with subsampled datasets in data/subsampled/. For example, to subsampled the *S. pneumoniae* pan-genome, run the following
```bash
cd scripts
./run_subsampling.sh spneumo_pangenome_k32
```
Furthermore, for generating negative queries, it is required to decompress [data/GRCh38.p14.chromosome1.prefix2M.fasta.xz](data/GRCh38.p14.chromosome1.prefix2M.fasta.xz) into experiments/01_build_and_query_memtime. Then run the experiment using
```bash
cd experiments/01_build_and_query_memtime
make
```
Notes:
- Since the resulting TSV tables are already in the repository, one needs to (re)move them to run the experiments.
- The number of cores provided to Snakemake can be changed in the Makefile (currently we use 4).
- The evaluated values of *k*, subsampling rates *r*, and datasets can all be changed in the [Snakefile](experiments/01_build_and_query_memtime/Snakefile).


## Figures

The [figures/](figures/) directory contains Fig. 2 and the associated supplementary plots, including the script used for their plotting.

## Contact

* Ondřej Sladký (ondra.sladky@gmail.com)
* Pavel Veselý (vesely@iuuk.mff.cuni.cz)
* Karel Břinda (karel.brinda@inria.fr)
