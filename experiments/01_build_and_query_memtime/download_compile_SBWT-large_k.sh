#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr SBWT_large_k/
	git clone https://github.com/algbio/SBWT.git SBWT_large_k
	(
		cd SBWT_large_k/
		git checkout b795178
		cd build/
		cmake .. -DCMAKE_CXX_COMPILER=g++-12 -DMAX_KMER_LENGTH=64
		make -j8
	)
)
ln -fs SBWT_large_k/build/bin/sbwt sbwt_large_k

./sbwt_large_k --help

echo "finished"
