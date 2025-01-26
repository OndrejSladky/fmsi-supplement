#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr ProphAsm/
	git clone git@github.com:prophyle/prophasm.git ProphAsm
	(
		cd ProphAsm/
		git checkout 017655e
		make -j 5
	)
)
ln -fs ProphAsm/prophasm 

./prophasm

echo "finished"
