#!/usr/bin/env bash
# This script should be run from the project root directory and will create the
# directory structure required for running the simulations pipeline
set -o pipefail -euv

# count number of sequences in fasta/fastq
cs () {
    if [ -f "$1" ]
    then
        case $1 in
            *.fastq)     wc -l "$1" | awk '{print $1/4}' ;;
            *.fq)        wc -l "$1" | awk '{print $1/4}' ;;
            *.fastq.gz)  zcat "$1" | wc -l | awk '{print $1/4}' ;;
            *.fq.gz)     zcat "$1" | wc -l | awk '{print $1/4}' ;;
            *.fasta)     grep -c "^>" "$1" ;;
            *.fa)        grep -c "^>" "$1" ;;
            *.fasta.gz)  zcat "$1" | grep -c "^>" ;;
            *.fa.gz)     zcat "$1" | grep -c "^>" ;;
            *)           echo "don't know how to count sequences in '$1'..." ;;
        esac
    else
        echo "'$1' is not a valid file"
    fi
}

data_dir="data/"
mkdir -p "$data_dir"

# Download the panx data
panx_dir="${data_dir}/all_gene_alignments"
tarball="${panx_dir}.tar.gz"
url="http://pangenome.tuebingen.mpg.de/dataset/Escherichia_coli/all_gene_alignments.tar.gz"

wget "$url" -O "$tarball"

# Extract panx data
mkdir -p "$panx_dir"
tar -xzf "$tarball" -C "$panx_dir" --strip-components=1
rm "$tarball"

# Remove amino acid alignments
case "$(uname -s)" in
    Linux*)     find "$panx_dir" -type f -name '*_aa_*' -print0 | xargs --null -I @ rm @ ;;
    Darwin*)    find "$panx_dir" -type f -name '*_aa_*' -print0 | xargs -I @ rm @ ;;
    *)          echo "Only Mac and Linux are supported"; exit 1;
esac

# duplicate sequence in files with only 1 sequence (clustalo requires more than one
# sequence in a file)
while IFS= read -r -d '' fasta
do
    num_seqs=$(cs "$fasta")
    if [ "$num_seqs" -eq 1 ]
    then
        zcat "$fasta" "$fasta" | awk '{print (/^>/ && (++c==2) ? ">duplicate_"substr($0, 2) : $0)}' | gzip -c - > temp
        mv temp "$fasta"
    fi
done <   <(find $panx_dir -name '*.fa.gz' -print0)
