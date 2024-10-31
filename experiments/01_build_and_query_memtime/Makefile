.PHONY: all help clean cleanall format test test1 cleantest build_programs

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

h=$(shell hostname)
$(info Hostname $h)

all: #build_programs
	snakemake -p -j 4 --rerun-incomplete --scheduler greedy --keep-going --rerun-triggers mtime --resources sbwtlimit=1

# Run test
test: cleantest #build_programs
	snakemake -p -j 20 --rerun-incomplete --config testing=True

# Run test with 1 thread
test1: cleantest #build_programs
	snakemake -p -j 1 --rerun-incomplete --config testing=True

format:
	snakefmt Snakefile

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -frv *_index/*
	rm -frv *_memtime/$h*
	rm -frv *_memtime/$h*
	rm -frv 42_query_fastas/*
	rm -frv 80_size_stats/$h*
	rm -frv 99_results/*$h*

cleantest:
	rm -fr 99_results/*$h__test.tsv

cleanall: clean ## Clean all
	rm -frv 00_hw_info/$h*


