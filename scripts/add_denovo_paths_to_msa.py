from pathlib import Path
import gzip

gene_paths = list(
    Path(snakemake.input.denovo_dir).glob(f"{snakemake.wildcards.gene}.*.fa")
)
new_msa_path = Path(snakemake.output[0])
original_msa_path = Path(snakemake.input.msa)
if original_msa_path.suffix == ".gz":
    open = gzip.open

if not new_msa_path.parent.is_dir():
    new_msa_path.parent.mkdir(parents=True, exist_ok=True)

with new_msa_path.open("w") as fh_out:
    original_msa = open(original_msa_path, "rt")
    fh_out.write(original_msa.read())

    for p in gene_paths:
        read_counter = 1

        with p.open() as fasta:
            for line in fasta:
                if line.startswith(">"):
                    fh_out.write(line.rstrip() + p.stem + f"_path{read_counter}\n")
                    read_counter += 1
                else:
                    fh_out.write(line)
