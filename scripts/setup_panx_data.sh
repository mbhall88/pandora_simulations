#!/usr/bin/env bash
# This script shhould be run from the project root directory and will create the
# directory structure required for running the simulations pipeline
data_dir="data/"

mkdir -p "$data_dir"

# Download the panx data
panx_dir="${data_dir}/all_gene_alignments"
tarball="${panx_dir}.tar.gz"
url="http://pangenome.tuebingen.mpg.de/dataset/Escherichia_coli/all_gene_alignments.tar.gz"

wget "$url" -O "$tarball"

# Extract panx data
tar -xzf "$tarball"
rm "$tarball"

# Remove amino acid alignments
find "$panx_dir" -type f -name '*_aa_*' -print0 | xargs -I @ rm @

# File alignments into subdirectories so there isn't too many files in one directory
files_per_dir=4000
output_subdir_num=0
n=0

for fpath in "$panx_dir"/*.fa.gz
do
    if [ "$n" -eq 0 ]; then
        outdir="${panx_dir}/${output_subdir_num}/"
        mkdir -p "$outdir"
        ((output_subdir_num++))
    fi
    mv "$fpath" "$outdir"
    ((n++))
    [ "$n" -eq "$files_per_dir" ] && n=0
done
