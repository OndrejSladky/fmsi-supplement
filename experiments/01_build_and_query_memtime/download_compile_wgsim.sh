#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr Wgsim/
	git clone https://github.com/lh3/wgsim.git Wgsim
	(
		cd Wgsim/
		git checkout a12da33
		gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
	)
)
ln -fs Wgsim/wgsim

./wgsim

echo "finished"
