#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr SSHash/
	git clone --recursive https://github.com/jermp/sshash.git SSHash
	(
		cd SSHash/
		git checkout d90ad37
		mkdir build
		cd build
		cmake ..
		make -j 5
	)
)
ln -fs SSHash/build/sshash

./sshash

echo "finished"
