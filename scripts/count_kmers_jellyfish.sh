#! /usr/bin/env bash

set -e
set -o pipefail
set -u
#!/bin/bash

FILE="$1"
K="$2"
OUT="$3"
TMP_FILE="tmp_$(basename $FILE).jf"


# Validate input arguments
if [[ -z "$FILE" || -z "$K" || -z "$OUT" ]]; then
  echo "Usage: $0 <fasta file> <k> <output file>"
  exit 1
fi

jellyfish count -m $K -s 2G -L 1 -C -t 3 -o "$TMP_FILE" "$FILE"
jellyfish stats "$TMP_FILE" | grep "Distinct:" | tr -cd [:digit:] >"$OUT"
rm -r "$TMP_FILE"

