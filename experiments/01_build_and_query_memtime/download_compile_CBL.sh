#! /usr/bin/env bash

set -e
set -o pipefail
set -u

PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x
echo "starting"

(
	rm -vfr CBL/
	git clone --recursive https://github.com/imartayan/CBL.git CBL
	(
		cd CBL/
		git checkout 328bcc6
		echo "======== compiling CBL for k=11 =============="
		RUSTFLAGS="-C target-cpu=native" K=11 PREFIX_BITS=20 cargo +nightly build --release --examples --target-dir target.k_11

		k=13
		echo "======== compiling CBL for k=$k =============="
		RUSTFLAGS="-C target-cpu=native" K=$k PREFIX_BITS=24 cargo +nightly build --release --examples --target-dir target.k_$k

		for k in {15..45..2}
		do
			echo "======== compiling CBL for k=$k =============="
			RUSTFLAGS="-C target-cpu=native" K=$k PREFIX_BITS=28 cargo +nightly build --release --examples --target-dir target.k_$k
		done
	)
)

echo "finished"
