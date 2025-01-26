#! /usr/bin/env bash

# set -e
# set -o pipefail
# set -u

#PS4='\[\e[32m\][$(date "+%Y-%m-%d %H:%M:%S") L${LINENO}]\[\e[0m\] '; set -x

# params: refGenome, mutation rate, seq. length, num. sequences, std. dev.
# output: stdout

WGSIM="../../wgsim/wgsim"
REF_GENOME=$1
MUT_RATE=$2
SEQ_LEN=$3
NUM_SEQS=$4
STD_DEV=$5

NUM_LINES=$((2*NUM_SEQS))

TMP_OUTPUT="/tmp/queries-$(uuidgen).fq"
TMP_OUTPUT_FA="/tmp/queries-$(uuidgen).fa"

"$WGSIM" -1"$SEQ_LEN" -s"$STD_DEV" -d0 -S42 -e0 -r"$MUT_RATE" -R0 -N"$NUM_SEQS" $REF_GENOME $TMP_OUTPUT /dev/null >/dev/null

seqtk seq -A -C "$TMP_OUTPUT" >"$TMP_OUTPUT_FA"
cat "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" "$TMP_OUTPUT_FA" | head -n $NUM_LINES
rm -rf "$TMP_OUTPUT" $TMP_OUTPUT_FA
