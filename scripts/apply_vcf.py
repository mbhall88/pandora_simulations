"""This script takes a VCF file and a fasta file and applies all variants in
that VCF to the fasta and writes the new sequence(s) to file.
"""
import pysam
import sys
import argparse
from pathlib import Path
from contextlib import ExitStack


def cli():
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--vcf",
                        type=str,
                        required=True,
                        help="VCF file.")

    parser.add_argument("--fasta",
                        type=str,
                        required=True,
                        help="Fasta file to apply variants to.")

    parser.add_argument("-o", "--output",
                        type=str,
                        help="""Path to output fasta file as. 
                        Default is STDOUT (-).""",
                        default="-")

    args = parser.parse_args()

    if args.output == "-":
        args.output = sys.stdout
    else:
        p = Path(args.output)
        assert p.parent.is_dir()
        args.output = p.open('w')

    assert Path(args.vcf).is_file()
    assert Path(args.fasta).is_file()

    return args


def apply_vcf_records_to_sequence(sequence, records):
    new_seq = ''
    pointer = 0
    for record in records:
        start = record.pos - 1
        end = start + len(record.ref)

        if start < pointer:
            raise IndexError(f"Variant coordinates overlap for {record.contig}")

        alt_idx = record.samples['sample']['GT'][0] - 1
        alt = record.alts[alt_idx]

        if sequence[start:end] != record.ref:
            raise ValueError(f"Ref. allele does not match for {record.contig}")

        new_seq += sequence[pointer:start] + alt
        pointer = end
    new_seq += sequence[pointer:]
    return new_seq


def main(args):
    with ExitStack() as stack:
        vcf = stack.enter_context(pysam.VariantFile(args.vcf))
        fasta = stack.enter_context(pysam.FastaFile(args.fasta))

        if vcf.index is None:
            raise IndexError("VCF file must be indexed with tabix.")

        for gene in fasta.references:
            sequence = fasta.fetch(reference=gene)

            if vcf.is_valid_reference_name(gene):
                records = list(vcf.fetch(contig=gene))
                new_seq = apply_vcf_records_to_sequence(sequence, records)
                print(f">{gene}\n{new_seq}", file=args.output)
            else:  # no variants for gene
                print(f">{gene}\n{sequence}", file=args.output)


if __name__ == '__main__':
    main(cli())
