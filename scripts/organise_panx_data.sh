#!/usr/bin/env bash
input="$1"
output="$2"
files_per_dir="$3"

find "$input" -type f -name '*_aa_*' -print0 | xargs -I @ rm @
outnum=0
n=0
for f in data/all_gene_alignments/*.fa.gz
do
    if [ "$n" -eq 0 ]; then
        outdir="$input"/"$outnum"
        mkdir -p "$outdir"
        ((outnum++))
    fi
    mv "$f" "$outdir"
    ((n++))
    [ "$n" -eq "$files_per_dir" ] && n=0
done
touch "$output"
