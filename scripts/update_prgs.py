import sys

sys.stderr = open(snakemake.log[0], "w")
from multiprocessing import Pool
from pathlib import Path
import logging

opts = snakemake.params.get("params", [])
prg_name = (lambda wildcards, output: Path(output.prg).with_suffix(""),)
msa_dir = Path(snakemake.input[0])
outdir = Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True, parents=True)

if snakemake.threads == 0:
    logging.info(f"0 processes requested. Using all available...")
    processes = None
else:
    processes = snakemake.threads


class MakePrgError(Exception):
    pass

include_genes = set(snakemake.params.genes)

def update_prgs(msa: Path, output: Path, gene: str):
    logging.debug(f"Updating PRG for {gene}...")
    prg_name = str(output).split(".")[0]
    args = ["make_prg", "from_msa", *opts, "-S", gene, "-n", prg_name, str(msa)]

    process = subprocess.Popen(
        args,
        stderr=sys.stderr,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise MakePrgError(
            f"Failed to execute make_prg for {gene} - check log file for error"
        )
    logging.debug(f"Finished updating PRG for {gene}")


jobs = []
for msa in msa_dir.glob("*.fa"):
    gene = msa.name.split(".")[0]
    if gene not in include_genes:
        continue
    prg_file = outdir / f"{gene}.prg"
    jobs.append((msa, prg_file, gene))


with Pool(processes=processes) as pool:
    logging.info(f"Updating {len(jobs)} PRGs...")
    pool.starmap(update_prgs, jobs)

logging.info("All PRGss updated!")
