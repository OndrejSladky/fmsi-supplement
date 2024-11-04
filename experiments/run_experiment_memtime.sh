#!/bin/bash

# get params.
if [ $# -ne 3 ]; then
    echo "Usage: $0 <file> <k> <subsampling rate>"
    echo "<prefix> = prefix of the FASTA file with the dataset"
    echo "      files ../data/<prefix>.fa.xz and ../data/<prefix>-refGenome.fa should exist"
    echo "<k> = k-mer size (positive integer)"
    echo "<subsampling rate> = rate so subsample k-mers; use 1 for no subsampling"
    exit 1
fi

k="$2"
r="$3"
g="$1"
file=../data/${g}.fa.xz

# dirs. for files created by the exp.
mkdir -p _query_fastas
mkdir -p _cbl_index
mkdir -p _fmsi_index
mkdir -p _bwa_index
mkdir -p _sbwt_index
mkdir -p _sshash_index
mkdir -p ../data/subsampled/

#############################################
echo "====== run subsampling if needed ======"
subsample=`echo "$r < 1" | bc -l`
if [ $subsample -eq 1 ]; then
    if [ -f ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz ]; then
        echo "../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz already exists"
    else
        echo "running subsampling for g=$g, k=$k, r=$r"
        ../scripts/subsample_kmers/subsample_kmers.py -k $k -r $r $file \
            | pv -l \
            | xz -9 -T4 \
            > ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz
        echo "created ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz"
    fi
else
    if [ ! -e "../data/subsampled/${g}_subsampled_k${k}_r1.0.fa.xz" ]; then
        ln -s ../../data/${g}.fa.xz "../data/subsampled/${g}_subsampled_k${k}_r1.0.fa.xz" 2>/dev/null && echo "created ../data/subsampled/${g}_subsampled_k${k}_r1.0.fa.xz"
    fi
    r=1.0
fi

#############################################
echo "====== compute_queries_streaming ======"
../wgsim/wgsim -1300 -d0 -S42 -e0 -r0 -R0 -N10000 ../data/${g}-refGenome.fa ${g}-queries-str.fq /dev/null
seqtk seq -a ${g}-queries-str.fq >_query_fastas/${g}-queries-str.fa
rm ${g}-queries-str.fq
      
#############################################
echo "====== compute_queries_negative ======="
../scripts/get_queries.py -k ${k} -cap 1000000 -print_header True -e <(xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz) <(xzcat ../data/GRCh38.p14.chromosome1.prefix2M.fasta.xz) >_query_fastas/${g}.k_${k}-queries-neg.fa
   
#############################################
echo "====== compute_queries_positive ======"
../scripts/get_queries.py -k ${k} -cap 1000000 -print_header True <(xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz) >_query_fastas/${g}.r_${r}.k_${k}-queries.fa
        
#############################################
echo "====== run_cbl ======"
xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz >_cbl_index/${g}_subsampled_k${k}_r${r}.fa
../scripts/benchmark.py --log "cbl_memtime_${g}.r_${r}.S_none.k_${k}.d_na.t_index.log" \
    "../CBL/target.k_${k}/release/examples/cbl build -c -o _cbl_index/${g}.r_${r}.S_none.k_${k}.d_na.fa.cblIndex _cbl_index/${g}_subsampled_k${k}_r${r}.fa"
rm -f _cbl_index/${g}_subsampled_k${k}_r${r}.fa

#############################################
echo "====== run_kmer_camel global ======"

../scripts/benchmark.py --log "camel_memtime_${g}.r_${r}.S_global.k_${k}.d_na.t_superstring.log" "../kmercamel/kmercamel -c -k ${k}  -p <(xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz) -a global >_fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na-no-opt.fa; ../kmercamel/kmercamel optimize -c -k ${k} -p _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na-no-opt.fa -a ones >_fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na.fa; rm _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na-no-opt.fa"

#############################################
echo "====== run_prophasm ======"
../scripts/benchmark.py --log "prophasm_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_superstring.log" \
    "../prophasm/prophasm -k ${k} -i <(xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz) -o _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa"

#############################################
echo "====== run_kmer_camel local d=1 ======"
../scripts/benchmark.py --log "camel_memtime_${g}.r_${r}.S_local.k_${k}.d_1.t_superstring.log" "../kmercamel/kmercamel -c -k ${k} -d 1 -p <(xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz) -a local >_fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1-no-opt.fa; ../kmercamel/kmercamel optimize -c -k ${k} -p _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1-no-opt.fa -a ones >_fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1.fa; rm _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1-no-opt.fa"

#############################################
echo "====== build_sbwt ======"
mkdir -p sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/
xzcat ../data/subsampled/${g}_subsampled_k${k}_r${r}.fa.xz >_sbwt_index/${g}_subsampled_k${k}_r${r}.fa
### TIME EFFICIENT VARIANT
# no streaming
../scripts/benchmark.py --log "sbwt_memtime_${g}.r_${r}.S_none.k_${k}.d_0.t_index.log" \
    "../SBWT/build/bin/sbwt build -k ${k} -t 1 -m 60 -d sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/ -i _sbwt_index/${g}_subsampled_k${k}_r${r}.fa -o _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_0.fa.sbwt --add-reverse-complements --variant plain-matrix --no-streaming-support"
rm -rf sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/*
# with streaming
../scripts/benchmark.py --log "sbwt_memtime_${g}.r_${r}.S_none.k_${k}.d_0.t_indexStreaming.log" \
    "../SBWT/build/bin/sbwt build -k ${k} -t 1 -m 60 -d sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/ -i _sbwt_index/${g}_subsampled_k${k}_r${r}.fa -o _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_0.fa-Streaming.sbwt --add-reverse-complements --variant plain-matrix"
rm -rf sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/*
### MEMORY EFFICIENT VARIANT
# no streaming
../scripts/benchmark.py --log "sbwt_memtime_${g}.r_${r}.S_none.k_${k}.d_1.t_index.log" \
    "../SBWT/build/bin/sbwt build -k ${k} -t 1 -m 60 -d sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/ -i _sbwt_index/${g}_subsampled_k${k}_r${r}.fa -o _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_1.fa.sbwt --variant rrr-split --no-streaming-support"
rm -rf sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/*
# with streaming
../scripts/benchmark.py --log "sbwt_memtime_${g}.r_${r}.S_none.k_${k}.d_1.t_indexStreaing.log" \
    "../SBWT/build/bin/sbwt build -k ${k} -t 1 -m 60 -d sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/ -i _sbwt_index/${g}_subsampled_k${k}_r${r}.fa -o _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_1.fa-Streaming.sbwt --variant rrr-split"
rm -rf sbwt-tmp/sbwt_tmp_${g}_subsampled_k${k}_r${r}.fa/
rm -f _sbwt_index/${g}_subsampled_k${k}_r${r}.fa

#############################################
echo "====== run_cbl_query_pos ======"
../scripts/benchmark.py --log "cbl_query_memtime_${g}.r_${r}.S_none.k_${k}.d_na.t_Pos.log" \
    "../CBL/target.k_${k}/release/examples/cbl query _cbl_index/${g}.r_${r}.S_none.k_${k}.d_na.fa.cblIndex _query_fastas/${g}.r_${r}.k_${k}-queries.fa >/dev/null"

#############################################
echo "====== compute_queries_negative_wRCs ======"
python3 ../scripts/add_RCs_to_fasta.py _query_fastas/${g}.k_${k}-queries-neg.fa _query_fastas/${g}.k_${k}-queries-neg-wRCs.fa

#############################################
echo "====== run_cbl_query_str ======"
../scripts/benchmark.py --log "cbl_query_memtime_${g}.r_${r}.S_none.k_${k}.d_na.t_Str.log" \
    "../CBL/target.k_${k}/release/examples/cbl query _cbl_index/${g}.r_${r}.S_none.k_${k}.d_na.fa.cblIndex _query_fastas/${g}-queries-str.fa >/dev/null"

#############################################
echo "====== run_sshash ======"
seqtk seq _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa >_sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa
../scripts/benchmark.py --log "sshash_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_index.log" \
    "../sshash/build/sshash build -i _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa -k ${k} -m 13 -o _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.sshash -s `echo $RANDOM`"
rm _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa

#############################################
echo "====== run_bwa ======"
../scripts/benchmark.py --log "bwa_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_index.log" \
    "../bwa/bwa index _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa"

#############################################
echo "====== run_fmsi ======"
../scripts/benchmark.py --log "fmsi_memtime_${g}.r_${r}.S_local.k_${k}.d_1.t_index.log" \
    "../fmsi/fmsi index -p _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1.fa" 

#############################################
echo "====== compute_queries_positive_wRCs ======"
python3 ../scripts/add_RCs_to_fasta.py _query_fastas/${g}.r_${r}.k_${k}-queries.fa _query_fastas/${g}.r_${r}.k_${k}-queries-wRCs.fa

#############################################
echo "====== run_cbl_query_neg ======"
../scripts/benchmark.py --log "cbl_query_memtime_${g}.r_${r}.S_none.k_${k}.d_na.t_Neg.log" \
    "../CBL/target.k_${k}/release/examples/cbl query _cbl_index/${g}.r_${r}.S_none.k_${k}.d_na.fa.cblIndex _query_fastas/${g}.k_${k}-queries-neg.fa >/dev/null"

#############################################
echo "====== run_fmsi ======"
../scripts/benchmark.py --log "fmsi_memtime_${g}.r_${r}.S_global.k_${k}.d_na.t_index.log" \
    "../fmsi/fmsi index -p _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na.fa" 

#############################################
echo "====== compute_queries_streaming_wRCs ======"
python3 ../scripts/add_RCs_to_fasta.py _query_fastas/${g}-queries-str.fa _query_fastas/${g}-queries-str-wRCs.fa

#############################################
echo "====== run_sshash_query_str ======"
../scripts/benchmark.py --log "sshash_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Str.log" \
    "../sshash/build/sshash query -i _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.sshash -q _query_fastas/${g}-queries-str.fa"
        
#############################################
echo "====== run_sbwt_query_neg ======"
### TIME EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_0.t_Neg.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_0.fa.sbwt -q _query_fastas/${g}.k_${k}-queries-neg.fa -o /dev/null"
### MEMORY EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_1.t_Neg.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_1.fa.sbwt -q _query_fastas/${g}.k_${k}-queries-neg-wRCs.fa -o /dev/null"

#############################################
echo "====== run_fmsi_query_neg ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_local.k_${k}.d_1.t_Neg.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1.fa -q _query_fastas/${g}.k_${k}-queries-neg.fa -O -s >/dev/null" 

#############################################
echo "====== run_fmsi_query_pos ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_global.k_${k}.d_na.t_Pos.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na.fa -q _query_fastas/${g}.r_${r}.k_${k}-queries.fa -O -s >/dev/null" 

#############################################
echo "====== run_bwa_query_pos ======"
../scripts/benchmark.py --log "bwa_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Pos.log" \
    "../bwa/bwa fastmap -l ${k} -w 999999 _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa _query_fastas/${g}.r_${r}.k_${k}-queries.fa >/dev/null"

#############################################
echo "====== run_fmsi_query_str ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_global.k_${k}.d_na.t_Str.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na.fa -q _query_fastas/${g}-queries-str.fa -O >/dev/null" 

#############################################
echo "====== run_sbwt_query_pos ======"
### TIME EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_0.t_Pos.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_0.fa.sbwt -q _query_fastas/${g}.r_${r}.k_${k}-queries.fa -o /dev/null"
### MEMORY EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_1.t_Pos.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_1.fa.sbwt -q _query_fastas/${g}.r_${r}.k_${k}-queries-wRCs.fa -o /dev/null"

#############################################
echo "====== run_sshash_query_pos ======"
../scripts/benchmark.py --log "sshash_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Pos.log" \
    "../sshash/build/sshash query -i _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.sshash -q _query_fastas/${g}.r_${r}.k_${k}-queries.fa"

#############################################
echo "====== run_fmsi_query_str ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_local.k_${k}.d_1.t_Str.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1.fa -q _query_fastas/${g}-queries-str.fa -O >/dev/null" 

#############################################
echo "====== run_sbwt_query_str ======"
### TIME EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_0.t_Str.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_0.fa-Streaming.sbwt -q _query_fastas/${g}-queries-str.fa -o /dev/null"
### MEMORY EFFICIENT VARIANT
../scripts/benchmark.py --log "sbwt_query_memtime_${g}.r_${r}.S_none.k_${k}.d_1.t_Str.log" \
    "../SBWT/build/bin/sbwt search -i _sbwt_index/${g}.r_${r}.S_none.k_${k}.d_1.fa-Streaming.sbwt -q _query_fastas/${g}-queries-str-wRCs.fa -o /dev/null"

#############################################
echo "====== run_fmsi_query_pos ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_local.k_${k}.d_1.t_Pos.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_local.k_${k}.d_1.fa -q _query_fastas/${g}.r_${r}.k_${k}-queries.fa -O -s >/dev/null" 

#############################################
echo "====== run_bwa_query_str ======"
../scripts/benchmark.py --log "bwa_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Str.log" \
    "../bwa/bwa fastmap -l ${k} -w 999999 _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa _query_fastas/${g}-queries-str.fa >/dev/null"

#############################################
echo "====== run_bwa_query_neg ======"
../scripts/benchmark.py --log "bwa_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Neg.log" \
    "../bwa/bwa fastmap -l ${k} -w 999999 _bwa_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.fa _query_fastas/${g}.k_${k}-queries-neg.fa >/dev/null"

#############################################
echo "====== run_sshash_query_neg ======"
../scripts/benchmark.py --log "sshash_query_memtime_${g}.r_${r}.S_prophasm.k_${k}.d_na.t_Neg.log" \
    "../sshash/build/sshash query -i _sshash_index/${g}.r_${r}.S_prophasm.k_${k}.d_na.sshash -q _query_fastas/${g}.k_${k}-queries-neg.fa"

#############################################
echo "====== run_fmsi_query_neg ======"
../scripts/benchmark.py --log "fmsi_query_memtime_${g}.r_${r}.S_global.k_${k}.d_na.t_Neg.log" \
    "../fmsi/fmsi query -p _fmsi_index/${g}.r_${r}.S_global.k_${k}.d_na.fa -q _query_fastas/${g}.k_${k}-queries-neg.fa -O -s >/dev/null" 
