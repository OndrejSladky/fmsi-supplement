#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr GGCAT/
	git clone https://github.com/algbio/ggcat.git GGCAT
	(
		cd GGCAT/
		git checkout fd59433
		cargo install --path crates/cmdline/ --locked
	)
)
ln -fs $HOME/.cargo/bin/ggcat

./ggcat -V

echo "finished"
