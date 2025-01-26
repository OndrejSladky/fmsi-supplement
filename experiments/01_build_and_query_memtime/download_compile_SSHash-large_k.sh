#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr SSHash_large_k/
	git clone --recursive https://github.com/jermp/sshash.git SSHash_large_k
	(
		cd SSHash_large_k/
		git checkout d90ad37
		mkdir build
		cd build
		cmake .. -D SSHASH_USE_MAX_KMER_LENGTH_63=On
		make -j 5
	)
)
ln -fs SSHash_large_k/build/sshash sshash_large_k

./sshash_large_k

echo "finished"
