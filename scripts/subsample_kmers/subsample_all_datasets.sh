#! /bin/bash

set -e
set -o pipefail
set -u

DATASETS_DIR=$1

for f in $DATASETS_DIR*.fa; do
    file=$(basename $f)
    echo ""
    echo "========================================"
    echo "STARTING PROCESSING ${file%.fa} $DATASETS_DIR"
    ./run_subsampling.sh "${file%.fa}" "$DATASETS_DIR"
    echo "========================================"
    echo ""
done
