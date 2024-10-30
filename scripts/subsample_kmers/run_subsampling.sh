#! /bin/bash

set -e
set -o pipefail
set -u

export OUTPUT_DIR="../../data/subsampled/"
mkdir -p $OUTPUT_DIR

GENOME=$1

function subsample {
    k="$2"
    r="$3"
    g="$1"
    if [ -f ${OUTPUT_DIR}${g}_subsampled_k${k}_r${r}.fa.xz ]; then
        echo "${OUTPUT_DIR}${g}_subsampled_k${k}_r${r}.fa.xz already exists"
        return
    fi
	echo "running subsampling for g=$g, k=$k, r=$r"
	./subsample_kmers.py -k $k -r $r ../../data/${g}.fa.xz \
		| pv -l \
		| xz -9 -T1 \
		> ${OUTPUT_DIR}${g}_subsampled_k${k}_r${r}.fa.xz
    echo "created ${OUTPUT_DIR}${g}_subsampled_k${k}_r${r}.fa.xz"
}

export -f subsample

s=""
for kk in 15 23 31
do
    ## the first rate is effectively just one, randomly chosen k-mer
    for rr in "0.1"
    do
		s+="$kk\n$rr\n"
	done
    if [ ! -e "${OUTPUT_DIR}${GENOME}_subsampled_k${kk}_r1.0.fa.xz" ]; then
        ln -s ../../data/${GENOME}.fa.xz "${OUTPUT_DIR}${GENOME}_subsampled_k${kk}_r1.0.fa.xz" 2>/dev/null && echo "created ${OUTPUT_DIR}${GENOME}_subsampled_k${kk}_r1.0.fa.xz"
    fi
done

printf $s \
	| parallel --max-args=2 -j1 subsample $GENOME
