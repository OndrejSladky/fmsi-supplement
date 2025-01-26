#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr SBWT/
	git clone https://github.com/algbio/SBWT.git SBWT
	(
		cd SBWT/
		git checkout b795178
		cd build/
		cmake .. -DCMAKE_CXX_COMPILER=g++-12 -DMAX_KMER_LENGTH=32
		make -j8
	)
)
ln -fs SBWT/build/bin/sbwt

./sbwt --help

echo "finished"
